#!/usr/bin/env python3

import sys
import itertools

if len(sys.argv) != 2:
    print("""Usage: deleteDuplicatesFasta.py myseqs.fa
    delete sequences with the same id and sequence from a fasta file
    if seqs are different, adds numeric suffix
    """)
else:
    myfile = sys.argv[1]

def parse_fasta(handle):
    """
    parses fasta file
    :param handle: filehandle of a fasta file
    :return: a generator, pairs of header/sequence items. Sequence devoid of newline characters.
    """
    fasta_iter = (list(g) for _,g in itertools.groupby(handle, lambda l: l.startswith('>')))
    for header_group, sequence_group in zip(*[fasta_iter]*2):
        header_string = ''.join(header_group).lstrip('>').rstrip()
        header = ''.join(header_group).lstrip('>').rstrip().split()[0]
        seq = ''.join([line.rstrip() for line in sequence_group])
        yield header, seq

record_dict = {}
with open(sys.argv[1], 'r') as f:
    for name, seq in  parse_fasta(f):
        if record_dict.get(name) is None:
            record_dict[name] = [seq]
        elif record_dict.get(name) == seq:
            pass
        elif record_dict.get(name) != seq:
            record_dict[name].append(seq)

for header, seq_ls in record_dict.items():
    n = 0
    for seq in seq_ls:
        n = n + 1
        print(f'>{header}_{n}')
        print(seq)



