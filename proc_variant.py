#! /usr/bin/env python

import re
import sys

min_variant_num = 3
seq_len = 48 

def remove_indel_bases(bases):
    
    # remaining = bases.replace('$', '')
    remaining = bases
    proc = ""
    while len(remaining) > 0:
        match = re.search(r'([\+\-])(\d+)', remaining)
        if match is None:
            proc = proc + remaining
            remaining = ""
        else:
            del_num = int(match.group(2))
            del_pos = match.start()

            proc = proc + remaining[0:del_pos]
            remaining = remaining[(del_pos + len(str(del_num)) + del_num + 1):len(remaining)]

    remaining = proc
    proc = ""
    while len(remaining) > 0:
        match = re.search(r'\^\S', remaining)
        if match is None:
            proc = proc + remaining
            remaining = ""
        else:
            pos = match.start()
            proc = proc + remaining[0:pos]
            remaining = remaining[(pos + 2):len(remaining)]

    proc = proc.replace('$', '')
    return proc


input_file = sys.argv[1]

with open(input_file, 'r') as hin:
    for line in hin:
        F = line.rstrip('\n').split('\t')
    
        # first simple check
        if int(F[3]) - F[4].count('.') - F[4].count(',') - F[4].count('>') - F[4].count('<') < min_variant_num:
            continue

        F4 = remove_indel_bases(F[4])
        if len(F4) != len(F[5]):
            print >> sys.stderr, "Inconsistent base and quality number"
            sys.exit(1)

        pos_vector = F[6].split(',')
        if len(F4) != len(pos_vector):
            print >> sys.stderr, "Inconsistent base and pos number"
            sys.exit()

        depth = 0
        base2num = {'A': 0, 'C': 0, 'G': 0, 'T': 0, 'N': 0, 'a': 0, 'c': 0, 'g': 0, 't': 0, 'n': 0}
        base2pos = {'A': [], 'C': [], 'G': [], 'T': [], 'N': [], 'a': [], 'c': [], 'g': [], 't': [], 'n': []}
        for i in range(len(F4)):
            if F4[i] in ['>', '<', '*']: continue
            depth = depth + 1
            if F4[i] == '.':
                base2num[F[2].upper()] = base2num[F[2].upper()] + 1
                base2pos[F[2].upper()].append(pos_vector[i])
            elif F4[i] == ',':
                base2num[F[2].lower()] = base2num[F[2].lower()] + 1 
                base2pos[F[2].lower()].append(pos_vector[i])
            else:
                base2num[F4[i]] = base2num[F4[i]] + 1
                base2pos[F4[i]].append(pos_vector[i])

        if depth == 0: continue


        alt = ''
        mis_rate = 0
        for nuc in ['A', 'C', 'G', 'T']:
            if nuc == F[2].upper(): continue
            cur_rate = float(base2num[nuc] + base2num[nuc.lower()]) / depth
            if cur_rate > mis_rate:
                mis_rate = cur_rate
                alt = nuc

        if mis_rate < 0.05: continue
        
        if len([x for x in base2pos[alt] + base2pos[alt.lower()] if int(x) >= 10 and int(x) <= seq_len - 10]) == 0: continue

        depth_p = base2num['A'] + base2num['C'] + base2num['G'] + base2num['T']
        depth_n = base2num['a'] + base2num['c'] + base2num['g'] + base2num['t']
        var_p = base2num[alt]
        var_n = base2num[alt.lower()]
        var_info = str(depth_p) + ',' + str(var_p) + ';' + str(depth_n) + ',' + str(var_n)
        strand_ratio = base2num[alt] / (base2num[alt] + base2num[alt.lower()])

        if F[1] == "7579311":
            pass

        if float(var_p + var_n) / int(F[3]) < 0.05: continue

        print '\t'.join([F[0], F[1], F[1], F[2]]) +'\t' +  alt + '\t' + str(round(mis_rate, 4)) + '\t' + str(depth_p + depth_n) + '\t' + \
            str(var_p + var_n) + '\t' + str(round(strand_ratio, 4)) +  '\t' + var_info
        


