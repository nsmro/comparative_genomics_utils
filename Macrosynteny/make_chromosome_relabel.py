#!/usr/bin/env python3

import argparse
import itertools

parser = argparse.ArgumentParser(description = 'prints 3 fields: chromosome relabel, original chromosome name, chromosome length')
parser.add_argument('genome',
                    help = 'genome, fasta format',
                    type = str)
parser.add_argument('-p', '--prefix',
                    help = 'prefix to add to chrom name (no separator)',
                    type = str,
                    required = True)
parser.add_argument('-n', '--number_chr',
                    help = 'number of chromosomes',
                    type = int,
                    required = True)
args = parser.parse_args()


def parse_fasta(filehandle):
    fasta_iter = (list(g) for _,g in itertools.groupby(filehandle, lambda l: l.startswith('>')))
    for header_group, sequence_group in zip(*[fasta_iter]*2):
        seqid = ''.join(header_group).split(None)[0].lstrip('>').strip()
        seq = ''.join(sequence_group)
        yield seqid,seq


chrom_list = []
with open(args.genome, 'r') as f:
    for name, chrom in parse_fasta(f):
        chrom_list.append([name, len(chrom)])

#sort chromosomes by descending lengths
x = 1
for old_name, length in sorted(chrom_list, key = lambda x: x[1], reverse = True):
    new_name = f'{args.prefix}{x:02}'
    print(f'{new_name}\t{old_name}\t{length}')
    x += 1
    if x > args.number_chr:
        break


