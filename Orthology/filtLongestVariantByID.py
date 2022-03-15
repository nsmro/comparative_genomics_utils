#! /usr/bin/env python3

import sys
import os
import re
import itertools
import argparse
import collections

fasta_record = collections.namedtuple('fasta_record',
                                       ['protein_acc', 'sequence', 'id'])

def parse_fasta(handle, ids_map = None):
    """
    parses fasta file, when given a handle
    ids_map is a dictionary
    :handle: filehandle of a fasta file
    :ids_map: dictionary to map gene ids to protein ids. If ids_map is None, gene id parsed from the fasta
    :return: a generator, pairs of header/sequence items. Sequence devoid of newline characters.
    """
    fasta_iter = (list(g) for _,g in itertools.groupby(handle, lambda l: l.startswith('>')))
    for header_group, sequence_group in zip(*[fasta_iter]*2):
        header_string = ''.join(header_group).lstrip('>').rstrip()
        prot_id = ''.join(header_group).lstrip('>').rstrip().split()[0]
        seq = ''.join([line.rstrip() for line in sequence_group])
        if ids_map is None:
            if args.filter_method == 'augustus': gene_re = r'(.*)\.\w+[^$, ^\s]'
            if args.filter_method == 'ensembl': gene_re = r'gene\:([^\s]+)'
            gene_id = re.search(gene_re, header_string).group(1)
        else:
            gene_id = ids_map.get(prot_id)
        if gene_id != None:
            yield fasta_record(prot_id, seq, gene_id)


def gff_to_id_map(handle, k_prot_id, k_gene_id):
    """
    :handle: gff filehandle
    :k_prot_id: key of the protein_id
    :k_gene_id: key of the gene_id
    :returns: a dictionary mapping prot_ids to gene_ids
    """
    output_dict = {}
    protID_re = k_prot_id+'\s*[:="\']*([^,|;"\']+)[;|"\s\n]*'
    gID_re = k_gene_id+'\s*[:="\']*([^,|;"\']+)[;|"\s\n]*'
    for line in handle:
        line = line.rstrip()
        try:
            gene_id = re.search(gID_re, line).group(1)
            prot_id = re.search(protID_re, line).group(1)
            if output_dict.get(prot_id) is None:
                output_dict[prot_id] = gene_id
        except AttributeError:
            pass #skips line if any re.search is "None"
    return output_dict


def write_longest(faiter, output_handle):
    """
    filters a fasta file, keeps only longest gene for each unique gene_id
    :faiter: fasta iterator, with 3 values: protein
    :output: name of the output file
    :returns: writes longest isoform to output file
    """
    by_gene_id = lambda gene: gene.id
    for _, variants in itertools.groupby(sorted(faiter,
                                                key = by_gene_id),
                                         key = by_gene_id):
        variant_ls = list(variants)
        single_variant = len(variant_ls) == 1
        if single_variant is True:
            long_v = variant_ls[0]
        else:
            long_v = sorted(variant_ls,
                            key = lambda v: len(v.sequence),
                            reverse = True)[0]
        output_handle.write(f'>{long_v.protein_acc}\n{long_v.sequence}\n')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description = """
        Filters a fasta file to keep the longest isoform of each gene,
        based on gene id.
        available filtering methods:
            `gff`: use gff to map gene_id to protein_ids.
                   default gff keys are the most common NCBI ones.
                   script should be able to parse most gff and some gtfs
            `ensembl`: parses gene id from ENSEMBL pep format.
            `augustus`: parses gene id from ENSEMBL pep format.
        """)
    parser.add_argument('--fasta', help = 'Common top-level parameter',
                        required=False)
    subparsers = parser.add_subparsers(help = 'method to filter longest isoforms',
                                       dest = 'filter_method')
    parent_parser = argparse.ArgumentParser(add_help = False) #circumvent help conflicy
    parent_parser.add_argument('-f', '--fasta', help = 'input fasta file', required = True)
    parent_parser.add_argument('-o',
                               '--output',
                                help = 'name of output file',
                                type = str)
    parser_gff = subparsers.add_parser('gff',
                                        help = 'use gff to map gene_id to protein_ids\
                                                default gff keys are the most common NCBI ones\
                                                script should be able to parse most gff and some gtfs',
                                        parents = [parent_parser])      
    parser_gff.add_argument('-a',
                            '--annot',
                            help ='gff/gtf file path',
                            required = True)
    parser_gff.add_argument('-kp',
                            '--key_protein_id',
                            help ='gff key, where "key_protein_id=accession",\
                                   default is "protein_id", for most NCBI datasets',
                            type = str,
                            default = 'protein_id')
    parser_gff.add_argument('-kg',
                            '--key_gene_id',
                            help ='gff key, found as "key_gene_id=geneid",\
                                   default is "gene", for most NCBI datasets',
                            type = str,
                            default = 'gene')
    parser_ensembl = subparsers.add_parser('ensembl',
                                           help = 'parses gene id from ENSEMBL pep format)',
                                           parents = [parent_parser])
    parser_augustus = subparsers.add_parser('augustus',
                                            help = 'parses gene id from augustus variants',
                                            parents = [parent_parser])
    args = parser.parse_args()

    if len(sys.argv) == 1:
        print("""
            Filters a fasta file to keep the longest isoform of each gene,
            based on gene id.
            available filtering methods:
                `gff`: use gff to map gene_id to protein_ids.
                       default gff keys are the most common NCBI ones.
                       script should be able to parse most gff and some gtfs
                `ensembl`: parses gene id from ENSEMBL pep format.
                `augustus`: parses gene id from ENSEMBL pep format.

            Try "filtLongestVariantByID.py --help" for more information
                 """)
        sys.exit()
    elif args.output is None:
        filename, file_extension = os.path.splitext(args.fasta)
        output_name = f'{filename}.longest_isoforms'
    else: output_name = args.output



    if args.filter_method == 'gff':
        with open(args.annot, 'r') as f:
            ids_map_dict = gff_to_id_map(handle = f,
                                         k_prot_id = args.key_protein_id,
                                         k_gene_id = args.key_gene_id)
    else:
        ids_map_dict = None
    with open(args.fasta, 'r') as g, open(output_name, 'w') as h:
        write_longest(faiter = parse_fasta(g, ids_map_dict),
                      output_handle = h)

