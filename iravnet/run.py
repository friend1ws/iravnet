#! /usr/bin/env python

from __future__ import print_function
import sys, subprocess, os, gzip, pkg_resources 

from .check_bam import *
from .proc_mpileup import *
from .filter_irav import *
from .annotate import *
from .validate import *

def get_main(args):

    genome_id, is_grc = check_refgenome(args.input_bam)

    if args.target_file is not None:
        target_file = args.target_file
    elif genome_id == "hg19":
        if is_grc:
            target_file = pkg_resources.resource_filename("iravnet", "data/boundary_proc.hg19.grc.bed")
        else:
            target_file = pkg_resources.resource_filename("iravnet", "data/boundary_proc.hg19.bed")
    elif genome_id == "hg38":
        if is_grc:
            target_file = pkg_resources.resource_filename("iravnet", "data/boundary_proc.hg38.grc.bed")
        else:
            target_file = pkg_resources.resource_filename("iravnet", "data/boundary_proc.hg38.bed")
    else:
        print("Output returned in the bam file check is strange.", file = sys.stderr)
        sys.exit(1)

    seqlen2count = check_seqlen(args.input_bam)

    if len(seqlen2count) == 0:
        print("No read in the input BAM file.", file = sys.stderr)
        sys.exit(1)

    """    
    if len(seqlen2count) > 1:
        print("Various lengths of reads in the input BAM file.", file = sys.stderr)
        sys.exit(1)
    """

    seqlen = list(seqlen2count.keys())[0]

    mpileup_command = ["samtools", "mpileup", args.input_bam, "-f", args.reference, "-l", target_file, "-q", str(args.min_mapq), "-O", "-o", args.output_file + ".tmp1"]
    # print(' '.join(mpileup_command))
    subprocess.check_call(mpileup_command, stderr = subprocess.DEVNULL)# , stdout = hout, stderr = subprocess.STDOUT)

    proc_mpileup(args.output_file + ".tmp1", args.output_file + ".tmp2", seqlen, args.min_variant_num, args.min_variant_ratio) 

    ir_ac_command = ["intron_retention_utils", "allele_count", args.input_bam, args.output_file + ".tmp2", args.output_file + ".tmp3", "--reference", args.reference]
    ir_ac_command = ir_ac_command + ["--donor_size", "3,6", "--acceptor_size", "6,1"]
    if genome_id == "hg38": ir_ac_command = ir_ac_command + ["--genome_id", "hg38"]
    subprocess.check_call(ir_ac_command)

    filter_irav(args.output_file + ".tmp3", args.output_file + ".tmp2", args.output_file)

    if not args.debug:
        subprocess.check_call(["rm", "-rf", args.output_file + ".tmp1"])
        subprocess.check_call(["rm", "-rf", args.output_file + ".tmp2"])
        subprocess.check_call(["rm", "-rf", args.output_file + ".tmp3"])



def annotate_main(args):

    annotate_vcf(args.input_vcf, args.output_vcf, args.gnomad_exome, args.gnomad_genome, args.clinvar)



def validate_main(args):

    hout = open(args.output_file, 'w')

    header_end_flag = False
    info_start_flag = False
    info_end_flag = False

    prefix = args.prefix

    with open(args.input_vcf, 'r') as hin:
        for line in hin:

            line = line.rstrip('\n')

            if line.startswith("##INFO") and not info_start_flag: info_start_flag = True

            if not line.startswith("##INFO") and info_start_flag and not info_end_flag:
                print('##INFO=<ID=' + prefix + '_DP,Number=1,Type=Integer,Description="Total read number">', file = hout)
                print('##INFO=<ID=' + prefix + '_AF,Number=1,Type=Float,Description="Allele Frequency">', file = hout)
                print('##INFO=<ID=' + prefix + '_AD,Number=1,Type=Integer,Description="Variant read number">', file = hout)
                info_end_flag = True

            if line.startswith('#'):
                print(line, file = hout)
            else:
                header_end_flag = True

            if not header_end_flag: continue

            F = line.split('\t')

            if F[0].startswith("#"): continue
            tchr, tpos, tref, tvar = F[0], F[1], F[3], F[4]

            tregion = tchr + ':' + tpos + '-' + tpos
            mpileup_command = ["samtools", "mpileup", args.input_bam, "-f", args.reference, "-r", tregion, "-q", str(args.min_mapq), "-o", args.output_file + ".tmp1.pileup"]
            subprocess.check_call(mpileup_command, stderr = subprocess.DEVNULL)

            tdepth, tvariant_num = validate_pileup(args.output_file + ".tmp1.pileup", tchr, tpos, tref, tvar)
            tvariant_ratio = float(tvariant_num) / float(tdepth) if tdepth != 0 else 0
        
            print(line + ";" + prefix + "_DP=" + str(tdepth) + ";" + prefix + "_AF=" + str(round(tvariant_ratio, 4)) + \
                  ";" + prefix + "_AD=" + str(tvariant_num), file = hout)

    hout.close()
    
    subprocess.check_call(["rm", "-rf", args.output_file + ".tmp1.pileup"])
            

def filt_bam_main(args):

    import annot_utils

    """
    gene_list = []
    with open(args.input_vcf, 'r') as hin:
        for line in hin:
            if line.startswith('#'): continue
            F = line.rstrip('\n').split('\t')
            INFO_list = F[7].split(';')
            for elm in INFO_list:
                if elm.startswith("GENE="):
                    gene = elm.replace("GENE=", '')
                    gene_list.append(gene)
    """
    gene2irav_info = {}
    with pysam.VariantFile(args.input_vcf) as vcf_h:
        for record in vcf_h.fetch():
            vchr, vpos = record.contig, record.pos
            vgene = record.info["GENE"]   
            gene2irav_info[vgene] = vchr, vpos


    genome_id, is_grc = check_refgenome(args.input_bam)
    annot_utils.gene.make_gene_info(args.output_bam + ".tmp.refGene.bed.gz", "refseq", genome_id, is_grc, False)

    gene2chr = {}
    gene2start = {}
    gene2end = {}
    """
    with gzip.open(args.output_bam + ".tmp.refGene.bed.gz", 'rt') as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')
            if F[3] in gene_list:
                gene2chr[F[3]] = F[0]
                if F[3] not in gene2start or int(F[1]) < int(gene2start[F[3]]):
                    gene2start[F[3]] = F[1]
                if F[3] not in gene2end or int(F[2]) > int(gene2end[F[3]]):
                    gene2end[F[3]] = F[2]
    """

    with gzip.open(args.output_bam + ".tmp.refGene.bed.gz", 'rt') as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')
            if F[3] in gene2irav_info:
                vchr, vpos = gene2irav_info[F[3]]
                if vchr == F[0] and int(vpos) >= int(F[1]) and int(vpos) <= int(F[2]):
                    gene2chr[F[3]] = F[0]
                    if F[3] not in gene2start or int(F[1]) < int(gene2start[F[3]]):
                        gene2start[F[3]] = F[1]
                    if F[3] not in gene2end or int(F[2]) > int(gene2end[F[3]]):
                        gene2end[F[3]] = F[2]

    region_list = []
    for gene in gene2chr:
        region = gene2chr[gene] + ':' + gene2start[gene] + '-' + gene2end[gene]
        region_list.append(region)


    # initialize the file
    hout = open(args.output_bam + ".tmp.unsorted.sam", 'w')
    hout.close()

    hout = open(args.output_bam + ".tmp.unsorted.sam", 'a')
    for region in sorted(region_list):
        subprocess.check_call(["samtools", "view", args.input_bam, region], stdout = hout, stderr = subprocess.DEVNULL)
    hout.close()

    hout = open(args.output_bam + ".tmp.unsorted.rmdup.sam", 'w')
    subprocess.check_call(["sort", "-u", args.output_bam + ".tmp.unsorted.sam"], stdout = hout)


    hout = open(args.output_bam + ".tmp.unsorted2.sam", 'w')
    subprocess.check_call(["samtools", "view", "-H", args.input_bam], stdout = hout, stderr = subprocess.DEVNULL)
    hout.close()

    hout = open(args.output_bam + ".tmp.unsorted2.sam", 'a')
    subprocess.check_call(["cat", args.output_bam + ".tmp.unsorted.rmdup.sam"], stdout = hout)
    hout.close()

    hout = open(args.output_bam + ".tmp.unsorted2.bam", 'w')
    subprocess.check_call(["samtools", "view", "-hbS", args.output_bam + ".tmp.unsorted2.sam"], stdout = hout)
    hout.close()

    hout = open(args.output_bam, 'w')
    subprocess.check_call(["samtools", "sort", args.output_bam + ".tmp.unsorted2.bam"], stdout = hout)
    hout.close()

    subprocess.check_call(["samtools", "index", args.output_bam])

    
    subprocess.check_call(["rm", "-rf", args.output_bam + ".tmp.unsorted.sam"])
    subprocess.check_call(["rm", "-rf", args.output_bam + ".tmp.unsorted.rmdup.sam"])
    subprocess.check_call(["rm", "-rf", args.output_bam + ".tmp.unsorted2.sam"])
    subprocess.check_call(["rm", "-rf", args.output_bam + ".tmp.unsorted2.bam"])
    subprocess.check_call(["rm", "-rf", args.output_bam + ".tmp.refGene.bed.gz"])
    subprocess.check_call(["rm", "-rf", args.output_bam + ".tmp.refGene.bed.gz.tbi"]) 
    


