#! /usr/bin/env python3

import argparse
import csv
import subprocess
import collections
from pathlib import Path


parser = argparse.ArgumentParser(description = """reads from a list of emapper annotation files,
some synteny results, a single clusters file
returns a tidy table with named columns, following names
multi_sp: id of the multi species block
block_id: id of the block
coordinates: coordinates with the format scaffold:start..end
node: prefix of the file where clusters is found, dot delimited
species: name of the species
acc_ls: comma-delimited list of accessions which show conserved synteny
name_ls: comma-delimited list of eggnog names, empty field if None)
nog_ls: comma-delimited list of NOG ids, empty field if None
nog_name_ls: comma-delimited list of NOG names, empty field if None
""")
parser.add_argument("emapper_folder",
                    help = "folder with eggnogmapper annotations. Only files with 22 fields will be parsed.")
parser.add_argument("-s",
                    "--synteny_files",
                    help = "synt file, output of MicroSynteny tool or SYNPHONI. Node name is guessed from first dot-delimited name of the filename",
                    required = True,
                    nargs = "+")
parser.add_argument("-c",
                    "--clusters_file",
                    help = "multi_species blocks, output of MicroSynteny tool or SYNPHONI",
                    required = True)
parser.add_argument("-o",
                    "--output",
                    help = "name of the outputfile",
                    default = "annotated_synteny.tsv")
args = parser.parse_args()


def load_multi_sp(filepath):
    """
    load multi_species nto a dict keys are block_ids, values are multi_sp ids
    """
    output_dict = {}
    with open(filepath, 'r') as f:
        for line in f:
            multi_sp, *block_ls = line.rstrip().split()
            output_dict.update({block : multi_sp for block in block_ls})
    return output_dict


def count_fields(filepath):
    """
    Counts the number of fields (tab-delimited) in a file using awk.
    requires subprocess
    """
    num = subprocess.check_output(f'head -n1 {filepath} | awk  --field-separator="\t" "{{ print NF }}"', shell = True)
    num = int(num.decode("utf-8").rstrip())
    return num


protein_annot = collections.namedtuple("protein_annot", "name nog nog_name")


def load_annot(filepath):
    field_check = count_fields(filepath) == 22
    if field_check is False:
        raise ValueError(f"the file {filepath} doesn't contain 22 fields. check input")
    output = {}
    with open(filepath, "r") as f:
        myreader = csv.reader(f, delimiter = "\t")
        for line in myreader:
            acc = line[0]
            acc_n = line[5]
            acc_nog = line[18].split("@")[0]
            acc_nog_name = line[21]
            output[acc] = protein_annot(name = acc_n, nog = acc_nog, nog_name = acc_nog_name)
    return output


def get_synt_node(filepath):
    """
    generator of block_id/node pairs
    """
    mypath = Path(filepath)
    node = mypath.stem.split(".")[0]
    with open(filepath, "r") as f:
        for line in f:
            block_id, *_ = line.rstrip().split()
            yield block_id, node


def process_synt(input_file, annot_map, node_map, multi_sp_map, processed_blocks, output_file):
    """process_synt writes info aboutthe file
    :param input_file: synt file name
    :param annot_map: keys are accessions, values are annotation_info
    :param node_map: keys are accessions, valuesare node_lists
    :param multi_sp_map: maps block ids to multi_species ids
    :param processed_blocks: blocks already written to output (set), edited inplace
    :param output_file: file where the stable will be saved
    """
    output_path = Path(output_file)
    with open(input_file, "r") as syntfile, open(output_path, "a") as annotation_table:
        syntreader = csv.reader(syntfile, delimiter = "\t")
        annot_writer = csv.writer(annotation_table, delimiter = "\t")
        for line in syntreader:
            block_id = line[0]
            if block_id not in processed_blocks:
                species = line[1]
                coords = line[7]
                acc_ls = line[9]
                node = ",".join(node_map.get(block_id))
                multi_sp = multi_sp_map[block_id]
                name_ls = []
                nog_ls = []
                nog_name_ls = []
                for acc in acc_ls.split(","):
                    acc_annot = annot_map.get(acc)
                    if acc_annot is None:
                        name_ls.append("")
                        nog_ls.append("")
                        nog_name_ls.append("")
                    else:
                        name_ls.append(acc_annot.name)
                        nog_ls.append(acc_annot.nog)
                        nog_name_ls.append(acc_annot.nog_name)
                row_ls = [multi_sp,
                          block_id,
                          coords,
                          node,
                          species,
                          acc_ls,
                          ",".join(name_ls),
                          ",".join(nog_ls),
                          ",".join(nog_name_ls)]
                annot_writer.writerow(row_ls)
                processed_blocks.add(block_id)

annot_folder = Path(args.emapper_folder)
file_ls = annot_folder.glob("*")
annot = {}
for annot_file in file_ls:
    annot |= load_annot(annot_file)        

multi_sp_dict = load_multi_sp(args.clusters_file)

synteny_files = args.synteny_files

node_blocks = collections.defaultdict(list)
for synt in synteny_files: # you can list as many input dicts as you want here
    for block_id, node in get_synt_node(synt):
        node_blocks[block_id].append(node)

header = "multi_sp\tblock_id\tcoordinates\tnode\tspecies\tacc_ls\tname_ls\tnog_ls\tnog_name_ls\n"

output_path = Path("test.tsv")
with open(output_path, "w") as f:
    f.write(header)

written_blocks = set()
for synt in synteny_files:
    process_synt(input_file = synt,
                 annot_map = annot,
                 node_map = node_blocks,
                 multi_sp_map = multi_sp_dict,
                 output_file = output_path,
                 processed_blocks = written_blocks)
