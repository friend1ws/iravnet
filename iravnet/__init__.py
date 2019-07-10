#! /usr/bin/env python

from .parser import create_parser

def main():

    parser = create_parser()
    args = parser.parse_args()
    args.func(args)

