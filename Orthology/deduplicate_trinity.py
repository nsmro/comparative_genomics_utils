#! /usr/bin/env python3

import sys
import itertools


if len(sys.argv) != 2:
    print("""Usage: deduplicatetrinity.py mytrinityfasta.fa > output.fa
        deduplicates trinity based on cluster IDs.
    """)
    sys.exit()



myfasta = sys.argv[1]

def parse_fasta(handle):
    """
    parses fasta file, when given a handle
    ids_map is a dictionary
    :handle: filehandle of a fasta file
    :return: a generator, pairs of header/sequence items. Sequence devoid of newline characters.
    """
    fasta_iter = (list(g) for _,g in itertools.groupby(handle, lambda l: l.startswith('>')))
    for header_group, sequence_group in zip(*[fasta_iter]*2):
        name = ''.join(header_group).lstrip('>').rstrip()
        cluster = "_".join(name.split()[0].split('_')[0:2])
        seq = ''.join([line.rstrip() for line in sequence_group])
        yield name, cluster, seq


with open(myfasta, 'r') as f:
    fasta_ls = sorted(list(parse_fasta(f)), key = lambda x: x[1])
    for cluster, group in itertools.groupby(fasta_ls, key = lambda x: x[1]):
        longest = sorted(group, key = lambda x: len(x))[-1]
        print(f">{longest[0]}\n{longest[2]}")

