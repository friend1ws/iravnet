#! /usr/bin/env python

import argparse

from .run import *
from .version import __version__

def create_parser():

    parser = argparse.ArgumentParser(prog = "savnet")

    parser.add_argument("--version", action = "version", version = "%(prog)s " + __version__)

    parser.add_argument("input_bam", metavar = "input.bam", default = None, type = str,
                        help = "the path to input bam file")

    parser.add_argument("output_file", metavar = "output.txt", default = None, type = str, 
                        help = "the path to the output")

    parser.add_argument("reference", metavar = "reference.fa", default = None, type = str,
                        help = "the path to the reference genome sequence")

    parser.add_argument("--gnomad_exome", type = str, default = None,
                        help = "the path to gnomad exome tabixed VCF file")

    parser.add_argument("--gnomad_genome", type = str, default = None,
                        help = "the path to gnomad whole genome tabixed VCF file")

    parser.add_argument("--min_variant_ratio", type = float, default = 0.05,
                        help = "the minimum value of the variant allele frequency for the variants")

    parser.add_argument("--min_variant_num", type = int, default = 3,
                        help = "the minimum number of suuporting reads for the variants")

    parser.add_argument("--min_mapq", type = int, default = 20,
                        help = "the threshould of mapping quality to count on the mpileup process for detecting candidate variants")

    """
    parser.add_argument("--effect_size_thres", type = float, default = 3.0,
                        help = "the thresould of effect size estimator used for simple edge pruning (default: %(default)s")

    parser.add_argument("--debug", default = False, action = 'store_true', help = "keep intermediate files")
    """
    
    return parser
 

