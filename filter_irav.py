#! /usr/bin/env python

import sys
import pysam

input_file = sys.argv[1]

gnomad_exome_db = pysam.Tabixfile("gnomad.exomes.r2.1.1.sites.vcf.bgz")
gnomad_genome_db = pysam.Tabixfile("gnomad.genomes.r2.1.1.sites.vcf.bgz")
# exac_db = pysam.Tabixfile("https://storage.googleapis.com/gnomad-public/release/2.1.1/vcf/genomes/gnomad.genomes.r2.1.1.sites.vcf.bgz")

with open(input_file, 'r') as hin:
    header2ind = {}
    header = hin.readline().rstrip('\n').split('\t')
    for i, cname in enumerate(header):
        header2ind[cname] = i

    print '\t'.join(header) + '\t' + "gnomAD_exome" + '\t' + "gnomAD_genome"

    for line in hin:
        F = line.rstrip('\n').split('\t')
        sj_n = int(F[header2ind["Splice_Junction_Negative"]])
        sj_p = int(F[header2ind["Splice_Junction_Positive"]])
        ir_n = int(F[header2ind["Intron_Retention_Negative"]])
        ir_p = int(F[header2ind["Intron_Retention_Positive"]])

        if ir_p < 3: continue
        if float(ir_p) / (ir_p + ir_n) < 0.9: continue
        if float(ir_p) / (ir_p + sj_p) < 0.9: continue

        chr_mut = F[header2ind["Chr_Mut"]]
        start_mut = int(F[header2ind["Start_Mut"]])  
        end_mut = int(F[header2ind["End_Mut"]])
        ref_mut = F[header2ind["Ref_Mut"]]
        alt_mut = F[header2ind["Alt_Mut"]] 

        AF_exome = 0
        for record_line in gnomad_exome_db.fetch(chr_mut, start_mut - 3, end_mut + 3):
            record = record_line.split('\t')
            if record[0] != chr_mut: continue
            if record[1] != str(start_mut): continue
            if record[3] != ref_mut: continue
            if record[4] != alt_mut: continue

            infos = record[7].split(';')
            for info in infos:
                if info.startswith("AF="):
                    AF_exome = float(info.replace("AF=", ''))
                            
                       
        AF_genome = 0
        for record_line in gnomad_genome_db.fetch(chr_mut, start_mut - 3, end_mut + 3):
            record = record_line.split('\t')
            if record[0] != chr_mut: continue
            if record[1] != str(start_mut): continue
            if record[3] != ref_mut: continue
            if record[4] != alt_mut: continue
            
            infos = record[7].split(';')
            for info in infos:
                if info.startswith("AF="):
                    AF_genome = float(info.replace("AF=", ''))

        print '\t'.join(F) + '\t' + str(round(AF_exome, 4)) + '\t' + str(round(AF_genome, 4))

