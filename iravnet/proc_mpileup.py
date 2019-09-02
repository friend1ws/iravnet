#! /usr/bin/env python

import re, sys, math


def proc_mpileup(input_file, output_file, seqlen, min_variant_num = 3, min_variant_ratio = 0.05, min_edge_margin = 10, min_pos_range = 5):

    # input_file = sys.argv[1]
    hout = open(output_file, 'w') 

    print('\t'.join(["Chr", "Start", "End", "Ref", "Alt", "Original_Ref", "Variant_Ratio", "Depth", 
                     "Variant_Num", "Strand_Ratio", "Variant_Num_Info"]), file = hout)

    with open(input_file, 'r') as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')
        
            # first simple check
            total_variant_num = int(F[3]) - F[4].count('.') - F[4].count(',') - F[4].count('>') - F[4].count('<') + F[4].count('+') + F[4].count('-')
            if total_variant_num < min_variant_num: continue
            if total_variant_num / int(F[3]) < min_variant_ratio: continue

            """
            F4 = remove_indel_bases(F[4])
            if len(F4) != len(F[5]):
                print >> sys.stderr, "Inconsistent base and quality number"
                sys.exit(1)

            pos_vector = F[6].split(',')
            if len(F4) != len(pos_vector):
                print >> sys.stderr, "Inconsistent base and pos number"
                sys.exit()
            """

            var2num = {}
            var2pos = {}
            var2num_plus = {}
            bases = F[4]
            # qualities = F[5].split('')
            positions = F[6].split(',')
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
                        var2pos[var].append(positions[base_ind])
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
                        # var_original = bases[0] + bases[2:(2 + indel_size)]
                        var_original = bases[0] + bases[(len(str(indel_size)) + 1):(len(str(indel_size)) + indel_size + 1)]
                        var = var_original.upper()
                        if var not in var2num:
                            var2num[var], var2pos[var], var2num_plus[var] = 0, [], 0
                        var2num[var] = var2num[var] + 1
                        var2pos[var].append(positions[base_ind])
                        if var == var_original: var2num_plus[var] = var2num_plus[var] + 1

                        # bases = bases[(2 + indel_size):]
                        bases = bases[(len(str(indel_size)) + indel_size + 1):]
                    base_ind = base_ind + 1

            if len(positions) != base_ind:
                print("Error???")
                sys.exit(1)

            if depth_p + depth_n == 0: continue


            bvar = ''
            bmis_rate = 0
            for var in var2num:
                if var2num[var] < min_variant_num: continue
                cur_rate = float(var2num[var]) / (depth_p + depth_n)
                if cur_rate > bmis_rate:
                    bmis_rate = cur_rate
                    bvar = var

            if bmis_rate < min_variant_ratio: continue
            if len([x for x in var2pos[bvar] if int(x) >= min_edge_margin and int(x) <= seqlen - min_edge_margin]) == 0: continue
            
            unique_positions = list(set([x for x in var2pos[bvar]]))
            if len(unique_positions) < min_variant_num: continue
            if max(unique_positions) - min(unique_positions) < min_pos_range: continue
            

            var_info = str(depth_p) + ',' + str(var2num_plus[bvar]) + ';' + str(depth_n) + ',' + str(var2num[bvar] - var2num_plus[bvar])
            strand_ratio = float(var2num_plus[bvar]) / var2num[bvar]

            start, end, ref, alt, original_ref = F[1], F[1], F[2], bvar, F[2]
            if bvar.startswith('-'): 
                ref = bvar[1:]
                alt = '-'
                start, end = str(int(start) + 1), str(int(start) + len(ref))
            if bvar.startswith('+'):
                ref = '-'
                alt = bvar[1:]

            print('\t'.join([F[0], start, end, ref, alt, original_ref, str(round(bmis_rate, 4)), str(depth_p + depth_n), str(var2num[bvar]),
                             str(round(strand_ratio, 4)), var_info]), file = hout)


    hout.close()


if __name__ == "__main__":

    import sys
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    seqlen = int(sys.argv[3])
    proc_mpileup(input_file, output_file, seqlen, min_variant_num = 3, min_variant_ratio = 0.05)
    
