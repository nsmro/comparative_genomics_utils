#1 /usr/bin/env python3

import ete3
import sys
import argparse


class clade:
    """
    A clade object. ingroups, sistergroup(s) and outgroups are isolated using the user-provided speciestree at instantiation.
    """
    def __init__(self, name, speciestree):
        """
        creates a clade object, from a species tree and a node name
        :param name: name of the noce where we define ingroups/outgroups
        :param speciestree: an ete3.Tree object, with nodenames
        :return: clade
        """
        self.name = name
        self.tree = speciestree
        self.ingroups = self.get_ingroups()
        self.outgroups = self.get_outgroups()
        self.novel_blocks = []
        self.ancestral_blocks = []
    def get_ingroups(self):
        """
        :return: nested list of ingroups of the clade
        """
        output_list = []
        for node in self.tree.traverse("preorder"):
            if node.name == self.name:
                for ingroup in node.children:
                    tmp_ig_list = []
                    for TaxonName in ingroup.iter_leaf_names(is_leaf_fn=None):
                        tmp_ig_list.append(TaxonName)
                    output_list.append(tmp_ig_list)
        return output_list
    def get_outgroups(self):
        """
        :return: flat list of outgroup species
        if the node is root, an empty list is returned
        """ 
        cached_tree = self.tree.copy()#create a copy in a different memory space of the species tree. node.detach() will affect cached_tree, but not the initial import of speciestree. For each new outgroup search, we'll bring back the copy of cached tree to its original state.
        for node in cached_tree.traverse("postorder"):
            if node.name == self.name:
                if node.is_root() is False:
                    node.detach()
                else:
                    return []
        return [name for name in cached_tree.iter_leaf_names(is_leaf_fn = None)]
    def block_content(self, species_map, threshold):
        """get_node_states edits in place a list of clade instances, adds multi species ids either to novel or ancestral
        :param multi_species_file: _description_
        :param list_clades: list of clade objects
        :param species_map: mapping of multi-species ids to species list
        :param threshold: species threshold. If block novel, needs to be found in 2 ingroups in at least n species
            if ancestral, block needs to be found in n species in ingroup, n species in outgroup
        :return: None, edits in place the clade objects
        """
        for multi_sp, species_ls in species_map.items():
            current_block_type = block_type(self, species_ls, threshold)
            if current_block_type == 'ancestral':
                self.ancestral_blocks.append(multi_sp)
            elif current_block_type == 'novel':
                self.novel_blocks.append(multi_sp)


def parse_block_info(block_file, multi_species_file):
    """parse_block_info 
    :param block_file: synt file
    :param multi_species_file: _description_
    :return: two dictionaries, one mapping the block_ids to  multi sp
    other mapping is mapping of cluster_id to species list who posess it
    """
    species_block = {}
    with open(block_file, "r") as f:
        species_block = {line.split()[0]: line.split()[1] for line in f}
    idx_map, species_map = {}, {}
    with open(multi_species_file, "r") as f:
        for line in f:
            multi_sp,  *idx_ls = line.rstrip().split()
            species_ls = [species_block[idx] for idx in idx_ls]
            species_map[multi_sp] = species_ls
            idx_map |= {idx: multi_sp for idx in idx_ls}
    return idx_map, species_map


def block_type(clade_object, species_ls, threshold):
    """
    Determines whether a block is ancestral or novel at a node of interest
    :param clade_object:  clade instance, node of interest
    :param species_ls: list of species
    :param threshold: species threshold. If block novel, needs to be found in 2 ingroups in at least n species
        if ancestral, block needs to be found in n species in ingroup, n species in outgroup
    :returns: either ancestral, novel or None
    """
    nb_species_per_ingroup = []
    for ingroup in clade_object.ingroups:
        nb_species_ingroup = sum(x in ingroup for x in species_ls)
        if len(ingroup) < threshold: #in the event that an ingroup is smaller than the specified threshold
            if nb_species_ingroup == len(ingroup): 
                nb_species_ingroup = threshold # if all the species of the ingroup posess the pair, threshold is satisfied
        nb_species_per_ingroup.append(nb_species_ingroup)
    nb_populated_ingroups = sum(x >= threshold for x in nb_species_per_ingroup)
    nb_species_outgroups = sum(x in clade_object.outgroups for x in species_ls)
    nb_species_ingroups = sum(nb_species_per_ingroup)
    block_is_novel = nb_populated_ingroups >= 2 and nb_species_outgroups == 0
    block_is_ancestral = any([nb_species_ingroups >= threshold and nb_species_outgroups >= threshold,
                         nb_species_ingroups >= threshold and nb_species_outgroups == len(clade_object.outgroups),
                         nb_populated_ingroups >= 2 and nb_species_outgroups != 0])
    if block_is_ancestral:
        return "ancestral"
    elif block_is_novel:
        return "novel"
    else:
        return None

    
def print_clusters_list(multi_sp_file, list_clade_instances, species_map, block_type):
    """print_blocks_list prints to stdout a list of blocks satisfying the specified paranmeters
    :param multi_sp_file: path of the multi species file file
    :param list_clade_instances: list of the nodes to analyze
    :param multi_sp_file: maps mutli_sp ids to list of species posessing the block
    :param block_type: ancestral/novel/total
    :returns: None, prints to STDOUT
    """
    with open(multi_sp_file,"r") as f:
        for line in f:
            multi_sp, *block_id_ls = line.rstrip().split('\t')
            block_id_str = "\t".join(block_id_ls) 
            node_list = []
            species_ls = set(species_map[multi_sp])
            for node in list_clade_instances:
                total_blocks = set(node.ancestral_blocks) & set(node.novel_blocks)
                conditions_to_keep_block = [
                    block_type == 'total' and multi_sp in total_blocks,
                    block_type == 'ancestral' and multi_sp in node.ancestral_blocks,
                    block_type == 'novel' and multi_sp in node.novel_blocks]
                if any(conditions_to_keep_block):
                    node_list.append(node.name)
            if node_list != []:
                node_ls_str = ','.join(set(node_list))
                print(f'{multi_sp}\t{node_ls_str}\t{",".join(species_ls)}\t{block_id_str}')


def print_blocks_list(blocks_file, list_clade_instances, idx_map, block_type):
    """print_blocks_list prints to stdout a list of blocks satisfying the specified paranmeters
    :param blocks_file: path of the synt file
    :param list_clade_instances: list of the nodes to analyze
    :param idx_map: maps block_ids to multi_sp ids
    :param block_type: ancestral/novel/total
    :returns: None, prints to STDOUT
    """
    with open(blocks_file, 'r') as f:
        for line in f:
            block_id, line_rest = line.rstrip().split('\t', maxsplit = 1)
            node_list = []
            multi_sp = idx_map[block_id]
            for taxon in list_clade_instances:
                total_blocks = set(taxon.ancestral_blocks) & set(taxon.novel_blocks)
                if block_type == 'total' and multi_sp in total_blocks:
                    node_list.append(taxon.name)
                if block_type == 'ancestral' and multi_sp in taxon.ancestral_blocks:
                    node_list.append(taxon.name)
                if block_type == 'novel' and multi_sp in taxon.novel_blocks:
                    node_list.append(taxon.name)
            if node_list != []:
                print(f'{multi_sp}\t{block_id}\t{line_rest}')


def print_short(list_clade_instances):
    """print_short report given a list of clade instances with node states
    :param list_clade_instances: list of clade objects, with block content inferred
    :returns: None, prints to stdout
    """
    ls_taxons = []
    for taxon in list_clade_instances:
        ls_taxons.append(taxon.name)
    print('\t'.join(['taxon', 'blocktype', 'count']))
    for taxon in list_clade_instances:
        print(f'{taxon.name}\tancestral\t{len(taxon.ancestral_blocks)}')
        print(f'{taxon.name}\tnovel\t{len(taxon.novel_blocks)}')


def print_tree(list_clade_instances, species_tree, block_type, report_type):
    """print_tree prints to STDOUT 
    :param list_clade_instances: list of clades where block content has been calculated
    :param species_tree: ete3.Tree object
    :param block_type: total, ancestral or novel
    :param report_type: _type_
    """
    cached_tree = species_tree.copy() 
    ls_taxons = [taxon.name for taxon in list_clade_instances]
    for node in cached_tree.traverse("preorder"):
        if node.name in ls_taxons:
            i = getattr(ls_taxons, 'index')(node.name) #ls_taxons and ls_clade_instances are essentially the same list order. We get the indes from ls_taxons, and use this index to extract class info from ls_clade_instances 
            print(i)
            node_clade = list_clade_instances[i]
            print(node_clade)
            if block_type == 'ancestral':
                node.name = f'{node.name}_{len(node_clade.ancestral_blocks)}'
            elif block_type == 'novel':
                node.name = f'{node.name}_{len(node_clade.novel_blocks)}'
            elif block_type == 'total':
                node.name = f'{node.name}_{len(node_clade.novel_blocks) + len(node_clade.ancestral_blocks)}'
    print(f"{block_type} block counts in {','.join(ls_taxons)}\n\n")
    if report_type == 'tree_ASCII':
        print(cached_tree.get_ascii(show_internal = True))
    elif report_type == 'tree_NH':
        print(cached_tree.write(format = 1))

parser = argparse.ArgumentParser(description = "Provides reports of block content of specified nodes")
parser.add_argument("-c",
                    "--multi_species_file",
                    help = "Tsv file, the multi-species blocks, output of microsynteny tool/SYNPHONI pipeline",
                    required = True)
parser.add_argument("-b",
                    "--blocks_file",
                    help = "Tsv file, microsyntenic blocks details, output of microsynteny tool/SYNPHONI pipeline",
                    required = True)
parser.add_argument("-s",
                    "--species_tree",
                    help = "Tree of the PREFIX in the second column of the block list,\
                        newick format with node names (e.g.\
                        (((D:0.723274,F:0.567784)E:0.067192,(B:0.279326,H:0.756049)B:0.807788); ).",
                    required = True)
parser.add_argument("-n",
                    "--node_names",
                    help = "Space separated list of names of nodes of interest from in the species tree\
                    (e.g. 'E B A C'.",
                    nargs = "+",
                    required = True)
parser.add_argument('-m',
                    '--species_threshold',
                    help = "Minimum number n of species per clade for node inference\
                    (e.g. for novel blocks, requires n species in at least two ingroups,\
                    and no species in outgroup\
                    for ancestral blocks, n species in one ingroup an n in one outgroup).\
                    By default, n=2. For ingroup/outgroup of size < n, all species of said\
                    ingroup/outgroup are required to posess the block\
                    (e.g. if ingroup size of 1, the species is required to posess.",
                    type = int,
                    default = 2)
parser.add_argument("-r",
                    "--report",
                    help = "Type of report printed to STDOUT [Default: short] : \n\
                            - 'short': number of blocks per node (default). \n\
                            - 'clusters_list': filters one multi-species block per line,\
                                cluster IDs (field 1), ancestral/novel nodes of the block (field 2),\
                                species list (field 3), block_ids (field 4+).\
                                For getting a *.clusters file with only a subset,\
                                pipe SyntByNode output to `cut -f1,4-`.\n\
                            - 'blocks_list': blocks within the filtered multi-species blocks\n\
                            - 'tree_ASCII': block count per specified node on an ASCII tree\n\
                            - 'tree_NH': block count per specified node on a newick tree",
                    choices = ["short", "clusters_list", "blocks_list", "tree_ASCII", "tree_NH"],
                    default = "short")
parser.add_argument("-t",
                    "--block_type",
                    help = "Specify whether you want to report all blocks ('total'),\
                        only the ones inherited from older nodes ('ancestral') or only 'novel' ones.",
                    choices = ["total","ancestral","novel"],
                    nargs = "+",
                    default = "total")
args = parser.parse_args()

if len(args.block_type) > 1 and args.report != "short":
    sys.stderr.write(f"multiple block types can be reported only when requesting 'short' report\
        for {args.report} require either ONLY total, ONLY ancestral or ONLY novel")
    sys.exit()


speciestree = ete3.Tree(args.species_tree, format=1)
node_ls = [clade(nodename, speciestree) for nodename in args.node_names]
idx_map, species_map = parse_block_info(args.blocks_file, args.multi_species_file)

for node in node_ls:
    node.block_content(species_map, 2)

if 'tree' in args.report:
    print_tree(args.node_names, speciestree)
elif args.report == 'short':
    print_short(args.node_names)
elif 'list' in args.report:
    sys.stderr.write(f'The {args.block_type} blocks are printed to stdout.\n')
    if args.report == 'clusters_list':
        print_clusters_list(args.multi_species_file, args.node_names, species_map)
    else:
        print_blocks_list(args.blocks_file, args.node_names, idx_map)

