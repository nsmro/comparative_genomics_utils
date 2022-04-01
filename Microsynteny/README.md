These are some scripts intended to deal with microsynteny files (Microsynteny tool, see Simakov et al. 2013).

For an example of a complete microsynteny run see [this page](https://github.com/nsmro/Syntenic_density_and_transitions/blob/main/Detailed%20methods/04.%20Microsynteny%20pipeline.md) 

# 1. List of the tools

## 1.1 bed2chrom.sh

Usage: `chrom2bed.sh bedfile PREFIX> chromfile`

Takes two positional arguments: a bed file (containing at least 6 fields UCSC standard), and a prefix (used to assign a species in the synteny pipeline)

The chromfile is printed to STDOUT

## 1.2. chrom2bed.sh

Usage: `chrom2bed.sh chromfile > bedfile`

Takes a single positional argument: a chrom file. Converts the chromfile to bed6 (UCSC standard) and prionts it to STDOUT

## 1.3 pymakeMap.py

Usage `pymakeMap.py -gff GFF_INPUT -p PREFIX -f FEATURE -k KEY -o OUTPUT_CHROM`

gff parser, converts to chrom format. Since chrom requires non redundant datasets, a fasta file or a clus file can be used to filter sequences to retain.

List of parameters:

* **-gff, --gff_input**: name of the input gff file

* **-p, --prefix**: prefix that'll be in the first column of the output chromfile

* **-f, --feature**: gff feature the lines of interest should posess in the third column of the gff (e.g. CDS)

* **-k**: key in the gff file where the gene accession that'll be printed in the ouput chrom file is loacated

* **-o, --output**: name of the output chromfile

* optional parameters:
  
  * **-d, --delete_redundancy**: transcripts of the gff are merged and considered isoforms of the same gene if there is exon overlap. Only the longest isoform of each gene is kept. If both transcripts have three exons or more, they get assigned to the same geneif they share one internal exon (not the first or last one). If one transcript has two exons (and the other by two or more), they get assigned to the same gene if they share two exon start/stops. Finally if one gene has only one exon, it gets merged if it overlaps with any other gene. All of this works in a overlap manner (if links are A--B and B--C, the group will comprise A,B,C even though A--C are not directly linked) . For optimal results, run the script in two steps. once on gff entries located on the "+" strand and another time on gff entries located on the "-" strand.
  
  * **-r, --ref_filter {fasta,clus}**: type of file used to filter the gff (keep only entries of interest with a **fasta** or a **clus** file)
  
  * **-F, --filter_file**: fasta or clus file used for filtering (use only with **-r**)
  
  * **-s, --suffix**: suffix found in the accessions of the fasta but absent of the gff, e.g. _1 (if you generated the proteins using transdecoder, for example). Use only if `-r fasta` and `-F` options are used

*example command*

`pymakeMap.py -gff HOMSA_example.gff -p HOMSA -f CDS -k Name -o HOMSA_example.chrom`

*example input gff*

```
NC_000001.11    Gnomon  CDS     399041  399100  .       -       0       ID=cds-XP_024307730.1;Parent=rna-XM_024451962.1;Dbxref=GeneID:112268260,Genbank:XP_024307730.1;Name=XP_024307730.1;gbkey=CDS;gene=LOC112268260;product=uncharacterized protein LOC112268260 isoform X1;protein_id=XP_024307730.1
NC_000001.11    Gnomon  CDS     379769  379870  .       -       0       ID=cds-XP_024307730.1;Parent=rna-XM_024451962.1;Dbxref=GeneID:112268260,Genbank:XP_024307730.1;Name=XP_024307730.1;gbkey=CDS;gene=LOC112268260;product=uncharacterized protein LOC112268260 isoform X1;protein_id=XP_024307730.1
NC_000001.11    Gnomon  CDS     373144  373323  .       -       0       ID=cds-XP_024307730.1;Parent=rna-XM_024451962.1;Dbxref=GeneID:112268260,Genbank:XP_024307730.1;Name=XP_024307730.1;gbkey=CDS;gene=LOC112268260;product=uncharacterized protein LOC112268260 isoform X1;protein_id=XP_024307730.1
NC_000001.11    Gnomon  CDS     365565  365692  .       -       0       ID=cds-XP_024307730.1;Parent=rna-XM_024451962.1;Dbxref=GeneID:112268260,Genbank:XP_024307730.1;Name=XP_024307730.1;gbkey=CDS;gene=LOC112268260;product=uncharacterized protein LOC112268260 isoform X1;protein_id=XP_024307730.1
NC_000001.11    Gnomon  CDS     358153  358183  .       -       1       ID=cds-XP_024307730.1;Parent=rna-XM_024451962.1;Dbxref=GeneID:112268260,Genbank:XP_024307730.1;Name=XP_024307730.1;gbkey=CDS;gene=LOC112268260;product=uncharacterized protein LOC112268260 isoform X1;protein_id=XP_024307730.1
```

*example returned line*

```
HOMSA    HOMSA_XP_02430773    NC_000001.11    -    3581    399100
```

## 1.4 pymakeRandChrom.py

This script is made to randomize a list of chrom files. Ideally, one would run the microsynteny pipleine one at least three sets of shuffled chrom files and see that in the shuffled genomes there's less than 5-10% of the number of blocks found in observed genomes (a finer grained orthology reduces noise).

What this script does is shuffle the accessions columns of a chrom file, and returns a shuffled version of the  chromfile, with an added suffix with rand and the number of the shuffle.

Input can be either a spece delimited list of files provided as posiitonal arguments, a newline-delimited list of chrom files piped to STDIN, or a file containing a newline delimited list of chrom files as positional argument. For the first two cases, the `-f` argument should be "chrom", for the third one, it should be "list".

The `-n` argument can be used to specify the number of independant chrom shuffling that should be performed (default is 1). the `-o`  argument is used to specify the ouput folder (if no folder is specified, shuffled chrom files are saved in the current working directory)

## 1.5 <u>BlocksByNode</u>.py

Both evolclust and MicroSynteny tool will return a single list of multi species blocks. Some blocks comprise only closely related species, so they're not really ancestral to most of the taxon sample.

The aim of this script is to filter multi-species blocks based on which taxonomic node they were present in (e.g. Was it found in the Last Common Ancestor of chordates? Is it also found outside of chordates). It requires the [ete python module >= 2.3](http://etetoolkit.org/docs/latest/tutorial/index.html) to parse a user-provided newick tree . Details about the procedure can be found in [this manuscript](https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-022-08304-2) (see Supplementary Fig. 2B to have a graphical summary)

The following arguments can be used:

* **-c, --clusters_id**: output of the MicroSynteny tool, tsv file. Contains mapping of block_ids to multi-specie blocks. I mostly used the ".clusters" extension.

* **-b, --block_list**: the synt file, output of the MicroSynteny tool, tsv file. 

* **-s, --species_tree**: Tree where the species PREFIX (as stated in the chrom files) are the leaves. Should be in newick format with node names, e.g`(((D:0.723274,F:0.567784)E:0.067192,(B:0.279326,H:0.756049)B:0.807788);`. While you *can* build such a tree from scratch, I'd recommend using [NCBI's Taxonomy browser to sketch the tree](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi) and then adjust it with [Mesquite](https://www.mesquiteproject.org/) according to the phylogenetic studies you trust the most. Polytomies are allowed.

* **-n, --node_names**: Space separated list of node names (as found in the species cladogram)

* **-m, --species_threshold**: Minimum number n of species per clade for determining if a block was present in a given node (for novel blocks, requires *m* species in at least two ingroups, and no species in outgroup. for ancestral/inherited blocks, requires *m* species in at least two ingroups, or *m* species in one ingroup and *m*  in the outgroup. By default, m = 2, considering that the monophyly of syntenic blocks isn't a reasonable assumption (see [this paper from Winter et *al*.](https://doi.org/10.1093/nar/gkw843)).

* **-r, --report**: type of report that should be printed to STDOUT. Multiple options can be specified.
  
  * **short**: number of blocks per node (default)
  
  * **clusters_list**: filtered multi-species blocks. Field 1 is the multi-species id, field 2 is the list of nodes specified in node names the block is found in (if the required block type is *novel*, then only one node will be in this field), field 3 and onward correspond to the input multi species file
  
  * **blocks_list**: blocks within the filtered multi-species blocks, adds a field to the left, the multi-species id
  
  * **tree_ascii**: adds the block count after the node name, and prints it to the terminal as an ASCII tree
  
  * **tree_ascii**: adds the block count after the node name, and prints it to the terminal as an newick (New Hampshire) tree

* **-t, --block_type**: type of block criteria used to filter:
  
  * **total**: all the blocks found in the specified node(s)
  
  * **novel**: only the novel blocks
  
  * **ancestral**: only the blocks inherited from an older node

## 1.6 make_synteny_table_emapper.py

This is to create a table summarizing the information about blocks, tidy format (one row per block).

The following arguments are required

* **emapper_folder (positional argument)**: folder where the annoptations by [eggnog-mapper](https://github.com/eggnogdb/eggnog-mapper/) are located.

* **-s, --synteny_files**: one or more synt file, output of the MicroSynteny tool. Multiple synt files can be used if the multi-species blocks of multiple synteny files are in the same **clusters** file

* **-c, --clusters_file**: multi-species file (clusters), output of the MicroSynteny tool

* **-o, --output**: name of the output file (default: annotated_synteny.tsv, in the current working directory)

The tidy table has the following fields:

* **multi_sp**: id of the multi-species block

* **block_id**: id of the block

* **coordinates**: coordinate string (chromosome:start..end)

* **node**: node parsed from the inputfile (first field of the filename, dot-delimited). 

* **species**: species prefix (as in the chrom file)

* **acc_ls**: comma-delimited list of accessions

* **name_ls**: comma-delimited list of predicted protein names (same order as accessions)

* **nog_ls**: comma-delimited list of NOG IDs (same order as accessions)

* **nog_name_ls**: comma-delimited list of NOG names (same order as accessions)

## 1.7 evolclust_prep.py

This is to prepare input files for [evolclust](https://github.com/Gabaldonlab/EvolClust). While it is possible to use the same orthology file as the Microsynteny tool, it's not optimal. While the MicroSynteny tool performs best with a finer grained orthology (e.g. the phylogenetic hierarchical orthogroups of [OrthoFinder](https://github.com/davidemms/OrthoFinder)), evolclust performs best with a coarser grained orthology (e.g. the homologs outputted by [HomoMCL](https://github.com/willpett/homomcl)).

This requires three arguments:

- **-c, --chrom_folder**: folder where all te chrom files are located. Only files withthe *chrom* extension will be globbed

- **-f, --family**: orthology file, clus format

- **-o, --output_prefix**: prefix of the output file

This returns three files. One **mcl** file (families, one line per family.) One **lst** file, where accessions are all renamed according to the evolclust format. One **pickle** file (HIGHEST_PROTOCOL), which contains a python dict (accessions of the **lst** file are mapped to original accessions)

## 1.8 evolclust_to_synt.py

This is to convert the [evolclust](https://github.com/Gabaldonlab/EvolClust) output to a synt/clusters file format, to filter blocks using scripts like Blocksbynode. Requires pickle file containing dict of aliases (created by evolclust_prep.py).

Argument are:

* **-c, --chrom_folder**: folder yhere the chrom files are located, will glob files woth chrom extension

* **-e,--evolclust_results**: the "cluster_families.complemented.txt" evolclust creates

* **-a,--alias_dict**: pickle file created by evolclust_prep.py

* **-o, --output**: name of the output files. One file wth synt extension and one file with clusters extension will be created

# 2. Description of internal formats

## 2.1 chrom format

**chrom** is a format that is used to define the location of genes on a given chromosome. chrom format has six tab-separated fields. 

The chrom fields are:

1. **prefix** - the species prefix of the animal (unique identifier to all the species of a given dataset, will appear in *synt ouput of the microsynteny pipelines)
2. **accession** - gene accession of the species.
3. **chrom** - The name of the chromosome (e.g. chr3, chrY, chr2_random) or scaffold (e.g. scaffold10671)
4. **strand** - Defines the strand. Either "." (=no strand), "+" or "-"
5. **start** - The starting position of the gene in the chromosome or scaffold. Can either be 1-based or 0-based. Only relative positions will be used to determine gene order
6. **end** - The ending position of the genein the chromosome or scaffold. Can either be 1-based or 0-based. Only relative positions will be used to determine gene order
* *example of a chrom file*: 

```
ADIVA    ADIVA_GSADVT00000001001    HG380758    -    410    2188
ADIVA    ADIVA_GSADVT00000003001    HG380758    +    7097    8229
ADIVA    ADIVA_GSADVT00000004001    HG380758    +    8892    9945
ADIVA    ADIVA_GSADVT00000005001    HG380758    +    10638    12990
```

## 2.2 clus format

**clus** is a format of orthogroups, with a varying number of fields. The first field is the orthogroup name, the second contains the number of proteins within the orthogroup, the third field and onwards contain accessions. See the Orthology folder to create clus files.

* *example of a clus file*

```
N0.HOG0000007    2    AURAU_Seg1367.3    EUPSC_cluster_18415
```

## 2.3 MicroSynteny output

Two files are generated by the Microsynteny tool. A synt file (information about individual blocks) and a multi-species block file (blocks grouped by their orthology).

## 2.3.1 Synt format

**synt** is a format that is used to define blocks with conserved synteny. synt format has at least 10 tab-separated fields.

The synt fields are:

1. **block_id** - the id of the current block (is used in the multi-species file to map block_ids to multi_species ids)
2. **species** - the species prefix of the animal (unique identifier in the first field of the chrom files)
3. **num_links** - number of blocks the current block is orthologous to
4. **linked_blocks** - comma-separated list of ids and species (e.g. `1234(SPECIES1),5678(SPECIES2)`) of the blocks the the current block is orthologous to
5. **num_linked_species** - number of species where orthologs to the current block are found
6. **linked_species** - comma-separated list of all the species where orthologs to the current block are found
7. **strand** - Defines the strand. Either "." (=no strand, all the genes are oriented differently), "+" (all the genes are on the plus strand) or "-" (all the genes are on the mnus strand)
8. **location** - coordinates of the blocks, as `chromosome:start..end`
9. **length** - length of the block in basepairs (i.e. end - start)
10. **gene_names** - comma-separated list of the accessions part of the syntenic block
11. **optional fields** - fields 11 and onward are not standard. They can be used for comments.

## 2.3.2 Multi-species format (clusters/blocks)

The multi species block format contains information about how the clustering of blocks into multi-species blocks (i.e orthologous blocks). This format has at least 2 tab-separated fields, but there's no maximal number of fields in a row.

The first field is the multi-species block-id (unique within a given multi-species file). The second field and onward comprise teh block_ids (as found in the associated synt file)
