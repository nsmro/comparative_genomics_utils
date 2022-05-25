#! /usr/bin/env python3

import argparse
import itertools
import collections
import ete3

parser = argparse.ArgumentParser(description = """
generates a 3 field dataframe with: multi_sp id, core ogs, and non-core ogs
Core orthogroups are orthogroups found in m species of 2+ ingroups or
m species ingroup and m species outgroup
""")
parser.add_argument('-s',
                    '--synt_file',
                    help = "name of the synt file",
                    required = True)
parser.add_argument('-mp',
                    '--multi_sp',
                    help = 'Multi species block (total), used to get multi_sp ID',
                    required = True)
parser.add_argument('-t',
                    '--species_tree',
                    help = 'tree used to infer ingroup and outgroup',
                    required = True)
parser.add_argument('-n',
                    '--node_name',
                    help = "node to which the multi_species blocks belong to",
                    required = True)
parser.add_argument('-m',
                    '--min_species',
                    help = 'min number of species to require a core og to be found (2+ ingroups or 1 ingroup and outgroup)',
                    default = 2,
                    type = int)
parser.add_argument('-c',
                    '--ortho',
                    help = 'orthology file, clus format',
                    required = True)
parser.add_argument('-o',
                    '--output',
                    help = 'output prefix', required = True)
args = parser.parse_args()


class clade:
    """
    A clade instance. As the class is created, we automatically isolate ingroups and outgroups using the user-provided speciestree.
    """
    def __init__(clade, name, speciestree):
        """
        creates a clade instance, from a species tree and a node name
        :param name: name of the noce where we define ingroups/outgroups
        :param speciestree: an ete3.Tree object, with nodenames
        """
        clade.name = name
        clade.ingroups = clade.get_ingroups(speciestree)
        clade.outgroups = clade.get_outgroups(speciestree)
        clade.ingroups_flat = list(itertools.chain.from_iterable(clade.ingroups))
    def get_ingroups(clade, speciestree):
        """
        :returns: ingroups of the clade in the user-provided species tree
        as a nested list
        """
        output_list = []
        for node in speciestree.traverse("preorder"):
            if node.name == clade.name:
                for ingroup in node.children:
                    tmp_ig_list = []
                    for TaxonName in ingroup.iter_leaf_names(is_leaf_fn=None):
                        tmp_ig_list.append(TaxonName)
                    output_list.append(tmp_ig_list)
        return output_list
    def get_outgroups(clade,speciestree):
        """
        :returns: outgroup as a flatlist
        if the node is root node empty list is returned
        """ 
        cached_tree = speciestree.copy()#create a copy in a different memory space of the species tree. node.detach() will affect cached_tree, but not the initial import of speciestree. For each new outgroup search, we'll bring back the copy of cached tree to its original state.
        for node in cached_tree.traverse("postorder"):
            if node.name == clade.name:
                if node.is_root() is False:
                    node.detach()
                else:
                    return []
        return [name for name in cached_tree.iter_leaf_names(is_leaf_fn = None)]


def load_ortho(filepath):
    """
    load orthology into a dict keys are accessions, values are ogs
    """
    output_dict = {}
    with open(filepath,'r') as f:
        for line in f:
            og, _, *acc_ls = line.rstrip().split()
            output_dict.update({acc:og for acc in acc_ls})
    return output_dict


def load_multi_sp(filepath):
    """
    load multi_speciesi nto a dict keys are block_ids, values are multi_sp ids
    """
    output_dict = {}
    with open(filepath, 'r') as f:
        for line in f:
            multi_sp, *block_ls = line.rstrip().split()
            output_dict.update({block : multi_sp for block in block_ls})
    return output_dict

species_tree = ete3.Tree(args.species_tree, format = 1)
myclade = clade(args.node_name, species_tree)

ortho_dict = load_ortho(args.ortho)
multi_sp_d = load_multi_sp(args.multi_sp)

#loads orthogroup content into block_dict
block_dict = {}
with open(args.synt_file, 'r') as f:
    for line in f:
        fields = line.rstrip().split()
        block_id = fields[0]
        species = fields[1]
        acc_ls = fields[9].split(',')
        og_ls = [ortho_dict[acc] for acc in acc_ls]
        multi_sp = multi_sp_d[block_id]
        if block_dict.get(multi_sp) == None:
            block_dict[multi_sp] = {}
            block_dict[multi_sp][species] = set(og_ls)
        elif block_dict[multi_sp].get(species) == None:
            block_dict[multi_sp][species] = set(og_ls)
        elif block_dict[multi_sp].get(species) != None:
            block_dict[multi_sp][species] = block_dict[multi_sp][species] | set(og_ls)


with open(args.output, "w") as f:
    for multi_sp, species_dict in block_dict.items():
        og_found_in = {} #list of species in which orthogroups are found
        total_og_set = set()
        for species, og_ls, in species_dict.items():
            for og in og_ls:
                total_og_set.add(og)
                if og_found_in.get(og) is None:
                    og_found_in[og] = set([species])
                else:
                    og_found_in[og].add(species)
        core_ogs = []
        for og, species_ls in og_found_in.items():
            og_ig_ls = [[sp for sp in ig if sp in species_ls] for ig in myclade.ingroups]
            nb_populated_ig = sum([1 if len(ig) >= args.min_species else 0 for ig in og_ig_ls ])
            nb_species_ig = sum([1 if sp in species_ls else 0 for sp in myclade.ingroups_flat])
            if nb_populated_ig >= 2:
                core_ogs.append(og)
            elif nb_species_ig >= args.min_species:
                populated_outgroup = sum([1 if sp in species_ls else 0 for sp in myclade.outgroups])
                if populated_outgroup >= args.min_species:
                    core_ogs.append(og)            
        non_core_ogs = total_og_set - set(core_ogs)
        f.write(f"{multi_sp}\t{','.join(core_ogs)}\t{','.join(non_core_ogs)}\n")



