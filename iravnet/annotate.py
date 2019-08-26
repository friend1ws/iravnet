#! /usr/bin/env python

import sys
import pysam

target_rnames = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X"]

def annotate_vcf(input_vcf, output_vcf, gnomad_exome, gnomad_genome):

    if gnomad_exome is not None: gnomad_exome_db = pysam.Tabixfile(gnomad_exome)
    if gnomad_genome is not None: gnomad_genome_db = pysam.Tabixfile(gnomad_genome)

    header_end_flag = False
    info_start_flag = False
    info_end_flag = False

    hout = open(output_vcf, 'w')
    with open(input_vcf, 'r') as hin:
        for line in hin:
            line = line.rstrip('\n')

            if line.startswith("##INFO") and not info_start_flag: info_start_flag = True

            if not line.startswith("##INFO") and info_start_flag and not info_end_flag:
                if gnomad_exome is not None: 
                    print('##INFO=<ID=GNOMAD_EXOME,Number=1,Type=Float,Description="The allele frequencies provided by gnomAD exome">', file = hout)
                if gnomad_genome is not None:
                    print('##INFO=<ID=GNOMAD_GENOME,Number=1,Type=Float,Description="The allele frequencies provided by gnomAD whole genome">', file = hout) 
                info_end_flag = True

            if line.startswith('#'):
                print(line, file = hout)
            else:
                header_end_flag = True

            if not header_end_flag: continue

            irav_record = line
            F = line.split('\t')

            if gnomad_exome is not None: 

                GNOMAD_EXOME = 0
                if F[0] in target_rnames:
                    for record_line in gnomad_exome_db.fetch(F[0], int(F[1]) - 3, int(F[1]) + 3):
                        record = record_line.split('\t')
                        if record[0] != F[0]: continue
                        if record[1] != F[1]: continue
                        if record[3] != F[3]: continue
                        if record[4] != F[4]: continue
                        infos = record[7].split(';')
                        for info in infos:
                            if info.startswith("AF="):
                                GNOMAD_EXOME = float(info.replace("AF=", ''))
                irav_record = irav_record + ";GNOMAD_EXOME=" + str(round(GNOMAD_EXOME, 4))     

            if gnomad_genome is not None:

                GNOMAD_GENOME = 0
                if F[0] in target_rnames:
                    for record_line in gnomad_genome_db.fetch(F[0], int(F[1]) - 3, int(F[1]) + 3):
                        record = record_line.split('\t')
                        if record[0] != F[0]: continue
                        if record[1] != F[1]: continue
                        if record[3] != F[3]: continue
                        if record[4] != F[4]: continue
                        infos = record[7].split(';')
                        for info in infos:
                            if info.startswith("AF="):
                                GNOMAD_GENOME = float(info.replace("AF=", ''))
                irav_record = irav_record + ";GNOMAD_GENOME=" + str(round(GNOMAD_GENOME, 4))

            print(irav_record, file = hout)

    hout.close()

