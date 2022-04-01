#! /usr/bin/env python3

import argparse
import pickle
from itertools import groupby
from pathlib import Path
from functools import reduce
import csv

parser = argparse.ArgumentParser(description = """
using a list of chrom files, and an orthology file and pickle file created by evolclöus_prep
Creates a synt and clusters file
""")
parser.add_argument("-c",
                    "--chrom_folder",
                    type = str,
                    required = True,
                    help = "folder where the chrom files are located. All the files with the *chrom extension will be extracted")
parser.add_argument("-e",
                    "--evolclust_results",
                    type = str,
                    required = True,
                    help = "the cluster_families.complemented.txt created by evolclust")
parser.add_argument("-a",
                    "--alias_dict",
                    type = str,
                    required = True,
                    help = "pickled file created by evolclust_prep.py alias dictionary")
parser.add_argument("-o",
                    "--output_prefix",
                    help = "prefix of the output files (auto mcl and lst extensions + _alias.tsv extension)",
                    default = "evolclust_resuts")
args = parser.parse_args()

def make_aliases(pickled_dict):
    """
    loads pickled dict of aliases 
    returns a reverse dict out of it
    """
    with open(pickled_dict, "rb") as f:
        names = pickle.load(f)
    return {v:k for k,v in names.items()}


def load_evol_clust(evolclust_results, aliases):
    """
    loads evolclust results, translates evolclust ids into original prots
    returns a nested dict: [multi_sp][block_id] -> protein lists
    """
    output = {}
    with open(evolclust_results, "r") as f:
        cluster_group = (list(g) for _,g in groupby(f, key = lambda l: l.startswith("#")))
        for header, blocks_str_group in zip(*[cluster_group]*2):
            multi_sp = header[0].lstrip("#").strip()
            output[multi_sp] = {}
            for blocks_str in blocks_str_group:
                block_id, prot_ls = blocks_str.split()
                block_id = "_".join((multi_sp, block_id))
                original_prots = [aliases[prot] for prot in prot_ls.split(";")]
                output[multi_sp][block_id] = original_prots
    return output


def load_prot_cords(chrom_file):
    """
    creates a nested dict
    [prot_id][feature] -> value
    feture = {"scaffold", "start", "end"}
    """
    output = {}
    with open(chrom_file, "r") as f:
        for line in f:
            _, accession, scaffold, _, start, end = line.rstrip().split('\t')
            output[accession] = {"scaffold": scaffold,
                                  "start": int(start),
                                  "end": int(end)}
    return output


def get_species(block_id):
    """
    helper func
    """
    return block_id.split("_")[2]


def write_clusters(filehandle, multi_species):
    """
    writes clusters based on synt
    """
    clusters_writer = csv.writer(filehandle, delimiter = "\t")
    for multi_sp, blocks_dict in multi_species.items():
        clusters_writer.writerow([multi_sp] + list(blocks_dict.keys()))


def write_synt(filehandle, multi_species, coords):
    """
    provided a filehandle, evolclust results and coords
    writes the blocks loaded into synt format
    """
    synt_writer = csv.writer(filehandle, delimiter = "\t")
    for multi_sp, blocks_dict in multi_species.items():
        num_links = len(blocks_dict.keys()) - 1
        block_species = list(set(map(get_species, blocks_dict.keys())))
        for block_id, acc_ls in blocks_dict.items():
            species = get_species(block_id)
            linked_blocks = ",".join((x for x in blocks_dict.keys() if x != block_id))
            linked_species = set(block_species)
            linked_species.remove(species)
            num_linked_species = len(linked_species)
            linked_species = ",".join(linked_species)
            strand = "."
            scaffold = set((coords[acc]["scaffold"] for acc in acc_ls)).pop()
            start = min((coords[acc]["start"] for acc in acc_ls))
            end = max((coords[acc]["end"] for acc in acc_ls))
            location = f"{scaffold}:{start}-{end}"
            length_block = end - start
            joined_acc = ",".join(acc_ls)
            synt_writer.writerow([block_id,
                                 species,
                                 num_links,
                                 linked_blocks,
                                 num_linked_species,
                                 linked_species,
                                 strand,
                                 location,
                                 length_block,
                                 joined_acc])

output_prefix = Path(args.output_prefix)
chromfile_folder = Path(args.chrom_folder)
chromfile_list = chromfile_folder.glob("*.chrom")

protein_aliases = make_aliases(args.alias_dict)
multi_sp_dict = load_evol_clust(args.evolclust_results, protein_aliases)
chrom_data = reduce(lambda a, b: {**a, **b},
                    map(load_prot_cords, chromfile_list))

with open(output_prefix.with_suffix(".synt"), "w") as f:
    write_synt(f, multi_sp_dict,chrom_data)

with open(output_prefix.with_suffix(".clusters"), "w") as f:
    write_clusters(f, multi_sp_dict)
