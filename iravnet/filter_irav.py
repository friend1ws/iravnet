#! /usr/bin/env python

import sys
import pysam

# input_file = sys.argv[1]

vcf_header = """\
##fileformat=VCFv4.3
##FILTER=<ID=PASS,Description="All filters passed">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total read number">
##INFO=<ID=AF,Number=1,Type=Float,Description="Allele Frequency">
##INFO=<ID=AD,Number=1,Type=Integer,Description="Variant read number">
##INFO=<ID=SB,Number=1,Type=Integer,Description="Strand bias (variant read number on the forward strand devided by total variant read number)">
##INFO=<ID=DPF,Number=1,Type=Integer,Description="Total read number on the forward strand">
##INFO=<ID=DPR,Number=1,Type=Integer,Description="Total read number on the reverse strand">
##INFO=<ID=ADF,Number=1,Type=Integer,Description="Variant read number on the forward strand">
##INFO=<ID=ADR,Number=1,Type=Integer,Description="Variant read number on the reverse strand">
##INFO=<ID=GENE,Number=1,Type=String,Description="Gene symbol of affected gene">
##INFO=<ID=MOTIF_POS,Number=1,Type=String,Description="Splicing motif detail">
##INFO=<ID=MOTIF_STRAND,Number=1,Type=String,Description="Strand of the affected gene">
##INFO=<ID=MOTIF_TYPE,Number=1,Type=String,Description="The type of affected splicing site (donor or acceptor)">
##INFO=<ID=SJ_WT,Number=1,Type=Integer,Description="The number of normally spliced reads with the reference allele">
##INFO=<ID=SJ_MT,Number=1,Type=Integer,Description="The number of normally spliced reads with the mutated allele">
##INFO=<ID=IR_WT,Number=1,Type=Integer,Description="The number of intron retained reads with the reference allele">
##INFO=<ID=IR_MT,Number=1,Type=Integer,Description="The number of intron retained reads with the mutated allele">
##INFO=<ID=IR_MT,Number=1,Type=Integer,Description="The number of intron retained reads with the mutated allele">\
"""

# gnomad_exome_db = pysam.Tabixfile("gnomad.exomes.r2.1.1.sites.vcf.bgz")
# gnomad_genome_db = pysam.Tabixfile("gnomad.genomes.r2.1.1.sites.vcf.bgz")
# exac_db = pysam.Tabixfile("https://storage.googleapis.com/gnomad-public/release/2.1.1/vcf/genomes/gnomad.genomes.r2.1.1.sites.vcf.bgz")

# def filter_irav(input_file, mut_file, output_file, gnomad_exome, gnomad_genome):
def filter_irav(input_file, mut_file, output_file):

    mutkey2info = {}
    with open(mut_file, 'r') as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')
            mutkey2info['\t'.join(F[:5])] = '\t'.join(F[5:])


    # if gnomad_exome is not None: gnomad_exome_db = pysam.Tabixfile(gnomad_exome)
    # if gnomad_genome is not None: gnomad_genome_db = pysam.Tabixfile(gnomad_genome)

    hout = open(output_file, 'w')
    print(vcf_header, file = hout)
    print('\t'.join(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"]), file = hout)
    
    with open(input_file, 'r') as hin:
        header2ind = {}
        header = hin.readline().rstrip('\n').split('\t')
        for i, cname in enumerate(header):
            header2ind[cname] = i

        # header_line = '\t'.join(header)
        # header_line = header_line + '\t' + '\t'.join(["Variant_Ratio", "Depth", "Variant_Num", "Strand_Ratio", "Variant_Num_Detail"])

        # if gnomad_exome is not None: header_line = header_line + '\t' + "gnomAD_exome"
        # if gnomad_genome is not None: header_line = header_line + '\t' + "gnomAD_genome"

        # print(header_line, file = hout)
        # print '\t'.join(header) + '\t' + "gnomAD_exome" + '\t' + "gnomAD_genome"

        for line in hin:
            F = line.rstrip('\n').split('\t')
            sj_n = int(F[header2ind["Splice_Junction_Negative"]])
            sj_p = int(F[header2ind["Splice_Junction_Positive"]])
            ir_n = int(F[header2ind["Intron_Retention_Negative"]])
            ir_p = int(F[header2ind["Intron_Retention_Positive"]])

            if ir_p < 3: continue
            if float(ir_p) / (ir_p + ir_n) < 0.9: continue
            # if float(ir_p) / (ir_p + sj_p) < 0.9: continue

            chr_mut = F[header2ind["Chr_Mut"]]
            start_mut = F[header2ind["Start_Mut"]]
            end_mut = F[header2ind["End_Mut"]]
            ref_mut = F[header2ind["Ref_Mut"]]
            alt_mut = F[header2ind["Alt_Mut"]] 

            motif_info = F[header2ind["Chr_Motif"]] + ':' + F[header2ind["Start_Motif"]] + '-' + F[header2ind["End_Motif"]]

            mut_key = '\t'.join([chr_mut, start_mut, end_mut, ref_mut, alt_mut])
            original_ref, variant_ratio, depth, variant_num, strand_ratio, variant_num_info = mutkey2info[mut_key].split('\t')
            variant_num_info_f, variant_num_info_r = variant_num_info.split(';')
            depth_f, variant_num_f = variant_num_info_f.split(',')
            depth_r, variant_num_r = variant_num_info_r.split(',')

            
            if alt_mut == "-":
                pos = str(int(start_mut) - 1)
                ref = original_ref + ref_mut
                alt = original_ref
            elif ref_mut == "-":
                pos = start_mut
                ref = original_ref
                alt = original_ref + alt_mut
            else:
                pos = start_mut
                ref, alt = ref_mut, alt_mut
            irav_line = '\t'.join([chr_mut, pos, '.', ref, alt, '.', '.'])

            irav_line = irav_line + '\t' + "DP=" + depth
            irav_line = irav_line + ';' + "AF=" + variant_ratio
            irav_line = irav_line + ';' + "AD=" + variant_num
            irav_line = irav_line + ';' + "SB=" + strand_ratio
            irav_line = irav_line + ';' + "DPF=" + depth_f
            irav_line = irav_line + ';' + "DPR=" + depth_r
            irav_line = irav_line + ';' + "ADF=" + variant_num_f
            irav_line = irav_line + ';' + "ADR=" + variant_num_r
            irav_line = irav_line + ';' + "GENE=" + F[header2ind["Gene_Symbol"]]
            irav_line = irav_line + ';' + "MOTIF_POS=" + motif_info
            irav_line = irav_line + ';' + "MOTIF_STRAND=" + F[header2ind["Strand_Motif"]]
            irav_line = irav_line + ';' + "MOTIF_TYPE=" + F[header2ind["Type_Motif"]] 
            irav_line = irav_line + ';' + "SJ_WT=" + str(sj_n)
            irav_line = irav_line + ';' + "SJ_MT=" + str(sj_p)
            irav_line = irav_line + ';' + "IR_WT=" + str(ir_n)
            irav_line = irav_line + ';' + "IR_MT=" + str(ir_p)


            print(irav_line, file = hout)


    hout.close()
    # if gnomad_exome is not None: gnomad_exome_db.close()
    # if gnomad_genome is not None: gnomad_genome_db.close()


