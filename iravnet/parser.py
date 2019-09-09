#! /usr/bin/env python

import argparse

from .run import *
from .version import __version__

def create_parser():

    parser = argparse.ArgumentParser(prog = "iravnet")

    parser.add_argument("--version", action = "version", version = "%(prog)s " + __version__)

    subparsers = parser.add_subparsers()

    ##########
    # get_mut
    get = subparsers.add_parser("get",
                                help = "Get intron retention associated variant")


    get.add_argument("input_bam", metavar = "input.bam", default = None, type = str,
                     help = "the path to input bam file")

    get.add_argument("output_file", metavar = "output.txt", default = None, type = str, 
                     help = "the path to the output")

    get.add_argument("reference", metavar = "reference.fa", default = None, type = str,
                     help = "the path to the reference genome sequence")

    get.add_argument("--target_file", type = str, default = None,
                     help = "the custom target file for exon-intron boundary")

    get.add_argument("--min_variant_ratio", type = float, default = 0.05,
                     help = "the minimum value of the variant allele frequency for the variants")

    get.add_argument("--min_variant_num", type = int, default = 3,
                     help = "the minimum number of suuporting reads for the variants")

    get.add_argument("--min_mapq", type = int, default = 20,
                     help = "the threshould of mapping quality to count on the mpileup process for detecting candidate variants")

    get.add_argument("--debug", default = False, action = 'store_true', help = "keep intermediate files")

    get.set_defaults(func = get_main)
    ##########

    ##########
    # annotation
    annotate = subparsers.add_parser("annotate",
                                     help = "Annotate variants obtained by get command")

    annotate.add_argument("input_vcf", metavar = "input.vcf", default = None, type = str,
                          help = "the path to the intron retention assocaiated variants list file obtained by get command")

    annotate.add_argument("output_vcf", metavar = "output.vcf", default = None, type = str,
                          help = "the path to the output")

    annotate.add_argument("--gnomad_exome", type = str, default = None,
                          help = "the path to gnomad exome tabixed VCF file")

    annotate.add_argument("--gnomad_genome", type = str, default = None,
                          help = "the path to gnomad whole genome tabixed VCF file")

    annotate.add_argument("--clinvar", type = str, default = None,
                          help = "the path to ClinVar tabixed VCF file")

    annotate.set_defaults(func = annotate_main)
    ##########

    ##########
    # validate 
    validate = subparsers.add_parser("validate",
                                     help = "Validate variants obtained by get command using a genomic BAM file")

    validate.add_argument("input_vcf", metavar = "input.vcf", default = None, type = str,
                          help = "the path to the intron retention assocaiated variants list file obtained by get comman")

    validate.add_argument("input_bam", metavar = "input.bam", default = None, type = str,
                          help = "the path to input bam file")

    validate.add_argument("output_file", metavar = "output.txt", default = None, type = str,
                          help = "the path to the output")

    validate.add_argument("reference", metavar = "reference.fa", default = None, type = str,
                          help = "the path to the reference genome sequence")

    validate.add_argument("--min_mapq", type = int, default = 20,
                          help = "the threshould of mapping quality to count on the mpileup process for detecting candidate variants")
 
    validate.add_argument("--prefix", type = str, default = "VAL",
                          help = "the prefix string added to the INFO keys in the output VCF file")

    validate.set_defaults(func = validate_main)
    ##########
    
    return parser
 

