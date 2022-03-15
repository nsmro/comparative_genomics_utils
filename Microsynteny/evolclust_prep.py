#! /usr/bin/env python3


import argparse
from operator import itemgetter
import itertools
from pathlib import Path
import pickle
import csv

parser = argparse.ArgumentParser(description = """
using a list of chrom files, and an orthology file
builds the necessary files for a evolclust run (list and mcl)
also creates alias file (original protein + id of evolclust run)
""")
parser.add_argument("-c",
                    "--chrom_folder",
                    type = str,
                    required = True,
                    help = "folder where the chrom files are located. All the files with the *chrom extension will be extracted")
parser.add_argument("-f",
                    "--family",
                    type = str,
                    required = True,
                    help = """this is a *.clus file (tsv file, one row per family, first field is family name, field 2 is 
                    nb of genes in family, fields 3 and onwards are accessions
                    """)
parser.add_argument("-o",
                    "--output_prefix",
                    reauired = True,
                    help = "prefix of the output files (auto mcl and lst extensions + _alias.pickle extension)")
args = parser.parse_args()

output_prefix = Path(args.output_prefix)
chrom_folder = Path(args.chrom_folder)
alias_dict = {}

for chromfile in chrom_folder.glob("*chrom"):
    current_chromfile = []
    chrom_id = 1
    chrom_protein = 1
    with open(chromfile, "r") as f:
        for line in f:
            prefix, acc, scaffold, _, start, _ = line.rstrip().split()
            start = int(start)
            current_chromfile.append([acc,scaffold, start])
    sorted_chromfile = sorted(current_chromfile, key = itemgetter(1))
    for chrom, group in itertools.groupby(sorted_chromfile, key = itemgetter(1)):
        sorted_chrom = sorted(group, key = itemgetter(2))
        acc_ls = [info[0] for info in sorted_chrom]
        acc_ls_alias = []
        for acc in acc_ls:
            idx = str(acc_ls.index(acc) + chrom_protein)
            acc_ls_alias.append(f"{prefix}_{chrom_id}_{idx}")
        chrom_id += 1
        chrom_protein += len(acc_ls)
        alias_dict.update(dict(zip(acc_ls, acc_ls_alias)))

#list file, positions in accessions
with open(output_prefix.with_suffix(".list"), "w") as f:
    for acc in alias_dict.values():
        f.write(acc + "\n")

#pickle dict of aliases
with open(output_prefix.with_suffix(".pickle"), "wb") as fp:
    pickle.dump(alias_dict, fp, protocol=pickle.HIGHEST_PROTOCOL)

#make new family with ids made for the evolclust run
with open(args.family, "r") as f, open(output_prefix.with_suffix(".mcl"), "w") as g:
    mywriter = csv.writer(g, delimiter = "\t")
    for line in f:
        line = line.rstrip().split()
        acc_ls = [alias_dict[acc] for acc in line[2:]]
        mywriter.writerow(acc_ls)

