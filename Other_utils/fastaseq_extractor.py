#! /usr/bin/env python3

import sys
import itertools
import argparse


parser = argparse.ArgumentParser(description='feed a sequence name, a fasta file, a start and a stop, returns portion of given sequence as fasta format')
parser.add_argument('-f', '--fasta',
					help = 'fasta file',
					required = True)
parser.add_argument('-i', '--sequence_id',
					help = 'full genome file',
					required = True)
parser.add_argument('-s', '--start',
					type = int,
					help = 'Use chrom for keeping only', default = '')
parser.add_argument('-e', '--end',
					type = int,
					help = 'length_file of the proteome')
args = parser.parse_args()


def parse_fasta(handle):
	"""
	fasta parse returns pairs of header/sequence items. Sequence devoid of newline characters.
	The fasta parser relies on the lazy evaluation of the generator.
	Zip takes one item in one generator, one item in the other (could change it in more),
	but since it's the same item that we call twice, we end up advancing the generator each time.
	"""
	fasta_iter = (list(g) for _,g in itertools.groupby(handle, lambda l: l.startswith('>')))
	for header_group, sequence_group in zip(*[fasta_iter]*2):
		header = ''.join(header_group).lstrip('>').rstrip().split()[0]
		seq = ''.join([part.strip(None) for part in sequence_group]) #remove trailing newlines
		yield header, seq

with open(args.fasta, 'r') as fastafile:
	for seqname, seq in parse_fasta(fastafile):
		if seqname == args.sequence_id:
			start = args.start -1
			end = args.end
			print(f'>{seqname}:{args.start}..{end}\n{seq[start:end]}')
			sys.exit()#exits once it finds the sequence of interest





