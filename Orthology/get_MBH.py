#! /usr/bin/env python3

import argparse
import Bio.Blast.Applications
import os
import pandas as pd
import itertools
import subprocess
import sys

parser = argparse.ArgumentParser(description = 'find mutual best hits from provided fasta files. If no anchor is provided, does all-against-all-pairwise comparisons. If results are already computed by a previous run, blasts won\'t be re-run')
parser.add_argument('fasta', help = 'Provide a sopace-separated list of fasta file as follows "PREFIX=path/to/myfile.fa". PREFIX being the prefix used for naming the output blast tables and mbh files', type = str, nargs = '+')
parser.add_argument('-a', '--anchor', help = 'anchor file e.g. when one wants to do blasts of a list of proteomes against a single anchor (e.g. for BLG labelling in macrosynteny plots)', type = str)
parser.add_argument('-r', '--blast_results', help = 'name of the folder where to save the blast results and *mbh hits', default = "MBH_results")
parser.add_argument('-t', '--num_threads', help = 'number of threads the blast should use', default = 1)
parser.add_argument('-op', '--only_prep', help = 'creates databases, only prints commands for slurm', action = 'store_true')
args = parser.parse_args()

if any(["=" in x for x in args.fasta]) is False:
    print("check how you provided files, respect the PREFIX=FILE.fasta syntax")
    sys.exit()


def do_fwd_rev_blast(prefix1, prefix2, dict_input):
    """
    :param dict_input: dict where keys are prefixes, values are paths to input proteomes
    :param prefix1: prefix of sp 1
    :param prefix2: prefix of sp 2
    :return: dict with keys are frozenset(prefix1,prefix2), values are location of blast results
    raw blast files written to disk in args.blast_results
    """
    output_dict = {}
    s1 = dict_input[prefix1]
    s2 = dict_input[prefix2]
    db1 = Bio.Blast.Applications.NcbimakeblastdbCommandline(dbtype = 'prot',
                                                            input_file = s1,
                                                            parse_seqids = True)
    db2 = Bio.Blast.Applications.NcbimakeblastdbCommandline(dbtype = 'prot', 
                                                            input_file = s2,
                                                            parse_seqids = True)
    if os.path.exists(f'{s1}.pdb') is False: #create db only if it's not already created
        pdb1 = subprocess.Popen(str(db1), shell = True)
        pdb1.wait()
    if os.path.exists(f'{s2}.pdb') is False:
        pdb2 = subprocess.Popen(str(db2), shell = True)
        pdb2.wait()
    fwd_out = f'{args.blast_results}/fwd.{prefix1}.{prefix2}.blastp.tblout'
    rev_out = f'{args.blast_results}/rev.{prefix1}.{prefix2}.blastp.tblout'
    output_dict[(prefix1, prefix2)] = {'fwd': fwd_out, 'rev': rev_out}
    fwd_blastp = Bio.Blast.Applications.NcbiblastpCommandline(query = s1,
                                                              db = s2,
                                                              out = fwd_out,
                                                              outfmt = '6',
                                                              max_target_seqs = 1,
                                                              num_threads = args.num_threads)
    rev_blastp = Bio.Blast.Applications.NcbiblastpCommandline(query = s2,
                                                              db = s1,
                                                              out = rev_out,
                                                              outfmt = '6',
                                                              max_target_seqs = 1,
                                                              num_threads = args.num_threads)
    if args.only_prep is True:
        print(fwd_blastp)
        print(rev_blastp)
    else:
        if os.path.exists(fwd_out) is False: #run blast only if it's not already done
            p1 = subprocess.Popen(str(fwd_blastp), shell = True)
            print(f'running forward blast on {prefix1}, {prefix2}')
            p1.wait()
        if os.path.exists(rev_out) is False:
            p2 = subprocess.Popen(str(rev_blastp), shell = True)
            print(f'running reverse blast on {prefix1}, {prefix2}')
            p2.wait()
    return output_dict


def make_mbh_file(prefix1, prefix2, tblout_blast_path_dict):
    """
    :param prefix1: prefix of sp 1
    :param prefix2: prefix of sp 2
    :param tblout_blast_path_dict: dict of blast tblout paths generated by do_fwd_rev_blast
    :return: nothing
    mbh files are written to disk in args.blast_results
    """
    fwd_filename = tblout_blast_path_dict[(prefix1, prefix2)]['fwd']
    rev_filename = tblout_blast_path_dict[(prefix1, prefix2)]['rev']
    fwd_results = pd.read_csv(fwd_filename, sep="\t", header = None)
    rev_results = pd.read_csv(rev_filename, sep="\t", header = None)
    headers = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore']
    fwd_results.columns = headers
    rev_results.columns = headers
    mbh_df = pd.merge(fwd_results, rev_results[['qseqid', 'sseqid']],
                left_on='sseqid', right_on='qseqid',
                how='outer')
    mbh_df = mbh_df.loc[mbh_df.qseqid_x == mbh_df.sseqid_y]# Discard rows that are not Mutual Best Hits
    mbh_df = mbh_df.groupby(['qseqid_x', 'sseqid_x']).max()# take maximum of each column
    output_file = f'{args.blast_results}/{prefix1}-{prefix2}.mbh'
    mbh_df[['qseqid_y','sseqid_y','bitscore']].to_csv(output_file, header = [prefix2, prefix1, 'bitscore'], sep='\t', index = False)
    print(f'{output_file} created')






input_proteomes_dict = {p.split('=', 1)[0]:p.split('=', 1)[1] for p in args.fasta}
prefix_ls = sorted(list(input_proteomes_dict.keys()))
if args.anchor is None:
    prefix_combinations = [(x,y) for x,y in itertools.combinations(prefix_ls, 2)]
else:
    anchor_prefix, anchor_file = args.anchor.split('=', 1)
    input_proteomes_dict.update({anchor_prefix:anchor_file})
    prefix_combinations = [(anchor_prefix, x) for x in prefix_ls]

if os.path.exists(args.blast_results) is False:
    os.makedirs(args.blast_results)


for prefix_a,prefix_b in prefix_combinations:
    dict_blast_results = do_fwd_rev_blast(prefix_a,prefix_b, input_proteomes_dict) #yes, the dict is updated at each iteration, but the values are only useful for one iteration
    if args.only_prep is False:
        make_mbh_file(prefix_a, prefix_b, dict_blast_results)

