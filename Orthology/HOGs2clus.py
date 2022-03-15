#! /usr/bin/env python3

import argparse

parser = argparse.ArgumentParser(description = 'Tranforms HOGs of orthofinder (v2.4) into clus format')
parser.add_argument('input', help = 'input file, found in the "/Phylogenetic_Hierarchical_Orthogroups" folder OrthoFinder version > 2.4')
parser.add_argument('output', help = 'output filename clus format')
args = parser.parse_args()


with open(args.input, 'r') as f, open(args.output, 'w') as g:
    _ = f.readline() # discard header line
    for line in f:
        HOG, _, _, *sp = line.strip().split()
        ouput_ls = [name.rstrip(',') for name in sp]# delete commas betyeen paralogs of the same HOG
        output_str = '\t'.join(ouput_ls)
        g.write(f'{HOG}\t{len(ouput_ls)}\t{output_str}\n')

