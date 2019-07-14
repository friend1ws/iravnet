#! /usr/bin/env python

import re, sys, math

def validate_pileup(input_file, tchr, tpos, tref, talt):

    with open(input_file, 'r') as hin:
        F = hin.readline().rstrip('\n').split('\t')

        if F[0] == '':
            return(0, 0)

        if F[0] != tchr or F[1] != tpos or F[2] != tref[0]:
            print("Error!!", file = sys.stderr)
            print(tchr, tpos, tref, talt)
            print(F, file = sys.stderr)
            sys.exit(1)

        var2num = {}
        var2num_plus = {}
        var2pos = {}
        bases = F[4]
        base_ind = 0
        depth_p, depth_n = 0, 0
        while bases != '':
            if bases[0] in ['>', '<', '*']: 
                base_ind = base_ind + 1
                bases = bases[1:]

            elif bases[0] in '^':
                bases = bases[2:]
            elif bases[0] in '$':
                bases = bases[1:]
            elif bases[0] in ['.', ',', 'A', 'C', 'G', 'T', 'N', 'a', 'c', 'g', 't', 'n']:
                if bases[0] not in ['.', ',']: 
                    var_original = bases[0]
                    var = var_original.upper()
                    if var not in var2num:
                        var2num[var], var2pos[var], var2num_plus[var] = 0, [], 0
                    var2num[var] = var2num[var] + 1
                    if var == var_original: 
                        var2num_plus[var] = var2num_plus[var] + 1

                if bases[0] in ['.', 'A', 'C', 'G', 'T', 'N']:
                    depth_p = depth_p + 1
                else:
                    depth_n = depth_n + 1

                bases = bases[1:]

                
                if len(bases) > 0 and bases[0] in ['+', '-']:

                    match = re.search(r'^[\+\-](\d+)', bases)
                    indel_size = int(match.group(1))
                    var_original = bases[0] + bases[2:(2 + indel_size)]
                    var = var_original.upper()
                    if var not in var2num:
                        var2num[var], var2pos[var], var2num_plus[var] = 0, [], 0
                    var2num[var] = var2num[var] + 1
                    if var == var_original: var2num_plus[var] = var2num_plus[var] + 1

                    bases = bases[(2 + indel_size):]

                base_ind = base_ind + 1


        # deletion
        if len(tref) > 1:
            bvar = '-' + tref[1:] 
        # insertion
        elif len(talt) > 1:
            bvar = '+' + talt[1:]
        # substitution
        else:
            bvar = talt

        vnum = var2num[bvar] if bvar in var2num else 0
        return(depth_p + depth_n, vnum)



if __name__ == "__main__":

    import sys
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    seqlen = int(sys.argv[3])
    proc_mpileup(input_file, output_file, seqlen, min_variant_num = 3, min_variant_ratio = 0.05)
    
