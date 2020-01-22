#! /usr/bin/env python

import sys
import pysam

target_rnames = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X"]

def annotate_vcf(input_vcf, output_vcf, gnomad_exome, gnomad_genome, clinvar):

    if gnomad_exome is not None: gnomad_exome_db = pysam.Tabixfile(gnomad_exome, encoding="utf-8")
    if gnomad_genome is not None: gnomad_genome_db = pysam.Tabixfile(gnomad_genome, encoding="utf-8")
    if clinvar is not None: clinvar_db = pysam.Tabixfile(clinvar, encoding="utf-8")

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
                if clinvar is not None:
                    print('##INFO=<ID=CLNV,Number=.,Type=String,Description="The Variation IDs of exactly matching ClinVar variants">', file = hout)
                    print('##INFO=<ID=CLNV_SIG,Number=.,Type=String,Description=Clinical significances for exactly matching variants">', file = hout)
                    print('##INFO=<ID=CLNV_MOTIF,Number=.,Type=String,Description="The Variation IDs of matching ClinVar variants at motif levels">', file = hout)
                    print('##INFO=<ID=CLNV_SIG_MOTIF,Number=.,Type=String,Description=Clinical significances for matching variants at motif levels">', file = hout)


            if line.startswith('#'):
                print(line, file = hout)
            else:
                header_end_flag = True

            if not header_end_flag: continue

            irav_record = line
            F = line.split('\t')

            if gnomad_exome is not None: 

                GNOMAD_EXOME = 0
                if F[0] in target_rnames + ["chr" + x for x in target_rnames]:
                    for record_line in gnomad_exome_db.fetch(F[0], int(F[1]) - 1, int(F[1]) + 1):
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
                if F[0] in target_rnames + ["chr" + x for x in target_rnames]:
                    for record_line in gnomad_genome_db.fetch(F[0], int(F[1]) - 1, int(F[1]) + 1):
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

            if clinvar is not None:
            
                CLNVs = []
                CLNV_SIGs = []
                if F[0] in target_rnames + ["chr" + x for x in target_rnames]:
                    for record_line in clinvar_db.fetch(F[0].replace("chr", ''), int(F[1]) - 1, int(F[1]) + 1):   
                        record = record_line.split('\t')
                        # if record[0] != F[0]: continue
                        if record[1] != F[1]: continue
                        if record[3] != F[3]: continue
                        if record[4] != F[4]: continue
                        CLNVs.append(record[2])
                        infos = record[7].split(';')
                        for info in infos:
                            if info.startswith("CLNSIG="):
                                CLNV_SIGs.append(info.replace("CLNSIG=", ''))

                CLNV = ','.join(CLNVs) if len(CLNVs) > 0 else "NA"
                CLNV_SIG = ','.join(CLNV_SIGs) if len(CLNV_SIGs) > 0 else "NA"
                irav_record = irav_record + ";CLNV=" + CLNV + ";CLNV_SIG=" + CLNV_SIG


                CLNV_MOTIFs = []
                CLNV_SIG_MOTIFs = []
                for info in F[7].split(';'):
                    if info.startswith("MOTIF_POS"):
                        mchr, mregion = info.replace("MOTIF_POS=", '').split(':')
                        mstart, mend = mregion.split('-')

                if mchr in target_rnames + ["chr" + x for x in target_rnames]:
                    for record_line in clinvar_db.fetch(mchr.replace("chr", ''), int(mstart) - 1, int(mend) + 1):
                        record = record_line.split('\t')
                        if int(record[1]) < int(mstart): continue
                        if int(record[1]) > int(mend): continue
                        CLNV_MOTIFs.append(record[2])
                        infos = record[7].split(';')
                        for info in infos:
                            if info.startswith("CLNSIG="):
                                CLNV_SIG_MOTIFs.append(info.replace("CLNSIG=", ''))

                CLNV_MOTIF = ','.join(CLNV_MOTIFs) if len(CLNV_MOTIFs) > 0 else "NA"
                CLNV_SIG_MOTIF = ','.join(CLNV_SIG_MOTIFs) if len(CLNV_SIG_MOTIFs) > 0 else "NA"
                irav_record = irav_record + ";CLNV_MOTIF=" + CLNV_MOTIF + ";CLNV_SIG_MOTIF=" + CLNV_SIG_MOTIF


                                 
            print(irav_record, file = hout)

    hout.close()

