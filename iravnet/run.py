#! /usr/bin/env python

from __future__ import print_function
import sys, subprocess, os, pkg_resources 

from .check_bam import *
from .proc_mpileup import *
from .filter_irav import *

def iravnet_main(args):

    genome_id, is_grc = check_refgenome(args.input_bam)

    if genome_id == "hg19":
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
    
    if len(seqlen2count) > 1:
        print("Various lengths of reads in the input BAM file.", file = sys.stderr)
        sys.exit(1)

    seqlen = list(seqlen2count.keys())[0]

    """
    hout = open(args.output_file + ".tmp1", 'w') 
    mpileup_command = ["samtools", "mpileup", args.input_bam, "-f", args.reference, "-l", target_file, "-q", str(args.min_mapq), "-O"]
    print(' '.join(mpileup_command))
    subprocess.check_call(mpileup_command, stdout = hout, stderr = subprocess.STDOUT)
    hout.close()

    proc_mpileup(args.output_file + ".tmp1", args.output_file + ".tmp2", seqlen, args.min_variant_num, args.min_variant_ratio) 

    ir_ac_command = ["intron_retention_utils", "allele_count", args.input_bam, args.output_file + ".tmp2", args.output_file + ".tmp3", "--reference", args.reference]
    if genome_id == "hg38": ir_ac_command = ir_ac_command + ["--genome_id", "hg38"]
    subprocess.check_call(ir_ac_command)
    """
    
    filter_irav(args.output_file + ".tmp3", args.output_file, args.gnomad_exome, args.gnomad_genome)
 
