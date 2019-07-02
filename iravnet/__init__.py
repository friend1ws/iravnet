#! /usr/bin/env python

from .parser import create_parser
from .run import iravnet_main

def main():

    parser = create_parser()
    args = parser.parse_args()
    iravnet_main(args)

