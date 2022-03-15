#! /usr/bin/python3

"""
usage: prepare_homomcl_length_file.py file.fasta
prints to stdout a single-line file with alternating fastaseq names and lengths in basepairs
"""

import sys
import itertools

fasta_file = sys.argv[1]


def parse_fasta(handle):
        """
        fasta parse returns pairs of header/sequence items. Sequence devoid of newline characters.
        The fasta parser relies on the lazy evaluation of the generator.
        Zip takes one item in one generator, one item in the other (could change it in more),
        but since it's the same item that we call twice, we end up advancing the generator each time.
        """
        fasta_iter = (list(g) for _,g in itertools.groupby(handle, lambda l: l.startswith('>')))
        for header_group, sequence_group in zip(*[fasta_iter]*2):
                header = ''.join(header_group).lstrip('>').strip()
                seq = ''.join([part.strip(None) for part in sequence_group]) #remove trailing newlines
                yield header, seq


with open(fasta_file, 'r') as f:
	for header, seq in parse_fasta(f):
		sys.stdout.write('{}\t{}\t'.format(header, len(seq)))
	sys.stdout.write('\n')
