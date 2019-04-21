#! /usr/bin/env python

import sys

input_file = sys.argv[1]

with open(input_file, 'r') as hin:
    header2ind = {}
    header = hin.readline().rstrip('\n').split('\t')
    for i, cname in enumerate(header):
        header2ind[cname] = i

    print '\t'.join(header)

    for line in hin:
        F = line.rstrip('\n').split('\t')
        sj_n = int(F[header2ind["Splice_Junction_Negative"]])
        sj_p = int(F[header2ind["Splice_Junction_Positive"]])
        ir_n = int(F[header2ind["Intron_Retention_Negative"]])
        ir_p = int(F[header2ind["Intron_Retention_Positive"]])

        if ir_p < 3: continue
        if float(ir_p) / (ir_p + ir_n) < 0.9: continue
        if float(ir_p) / (ir_p + sj_p) < 0.9: continue

        print '\t'.join(F)


