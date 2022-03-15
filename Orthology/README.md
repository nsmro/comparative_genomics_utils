These are some scripts for orthology/phylogeny. The format used for storing orthology information (clus format) is described in microsynteny section..

# 1. List of the tools

## 1.1 filtLongestVariantByID.py

This is a script to filter fasta files using a gff. Aim is to keep the longest isoform of each gene. Since the MicroSynteny tool treats isoforms as paralogous genes, it's best to filter isoforms before orthology inference and microsynteny detection.

The script has three available methods for filtering the input fasta file: `gff`, `ensembl` or `augustus`. In every case, only the longest isoform of each respective gene will be kept. If the dataset you are using does not include gene ids but is redundant, the  `-d` flag of `microsynteny/pyMakeMap.py` might be another option to filter transcripts.

Each mode can be called using the filtering method method as a positional argument: e.g `filtLongestVariantByID.py gff <mode-specific arguments>`

The arguments are as follows:

* **-f, --fasta**: fasta file to filter (parameter common to all the filering methods).

* **-o, --output**: name of the ouput fasta where only the longest isoforms are kept (parameter common to all the filering methods). 

* **gff** filtering method. This method should be used if gene ids are present in the gff, in the same line as the protein ids (i.e. name of the sequences in the fasta file input). Only the longest isoform of each gene will be kept:
  
  * **-a, --annot**: name of the gff file
  
  * **-kp, --key_protein_id**: id of the protein (i.e. found in the headers of the fasta input). Default is "protein_id", the most common found in the gff3 files from GeneBank.
  
  * **-kg, --key_gene_id**: id used to group the isoform of the same gene. Default is "gene_id", most common in gff3 files from GeneBank

* **ensembl** method.  This takes advantage of the fact that ensembl fasta files include gene ids in the fasta headers. No additional argument required.

* **augustus** method. Fasta names are treated as following this format: ">gene_id.isoform_id". First field (dot-delimited) is treated as gene id, second field (dot-delimited) is treated as isoform id. Only longest isoform of its respective gene is kept

## 1.2 deduplicate_trinity.py

This is a script to deduplicate transcriptomes assembled using Trinity.

Usage: `deduplicatetrinity.py mytrinityfasta.fa > output.fa`

Say our file comprises the following headers: `comp0_c0_seq1`, `comp0_c0_seq2` and `comp0_c0_seq3`

The sequences are considered as isoforms of the same gene based on the first two undescore-delimited fields (here: "comp0_c0"). Only the longest will be kept.

## 1.3. get_MBH.py

This is a wrapper script for multiple blastp calls. This requires the [BLAST+](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/) suite and the [Biopython](https://biopython.org/wiki/Download) and [pandas](https://pandas.pydata.org/getting_started.html) packages. Either performs an all-against-all pairwise blast of multiple proteomes, or blasts multiple files against an anchor proteome. This creates mbh files (mutual best hits, comprises three fields: accession1, accession2, and the bitscore).

In the case of an all-against-all blast with not too many organisms (no real threshold there, depends how closely they are related), all the mbh files can be concatenated, loaded into a graph analysis API (igraph in R, networkx in python) since an mbh file is essentially an edge list. Finallyone can look for cliques in this graph to approximate orthogroups (again, not too many species).

 In the case of a blast against an anchor (reference transcriptome), one could fuse mbh using the anchor as a key, and retain only the proteins that have a mutal best hit in all the species to approximate orthogroups.

The parameters of the script are as follows:

* **fasta (positional)**: a space-separated list of fasta files, formated such as: "PREFIX=path/to/my/file.fa". PREFIX is the string that will be used to name the output blast and mbh tables. 

* **-a, --anchor**: name of the anchor file, also formated as: "PREFIX=path/to/my/anchor.fa". PREFIX is the string that will be used to name the output blast and mbh tables

* **-r, --blast_results**: name of the folder for the mbh files. If the folder doesn't existm it's created

* **-t, --num_threads**: number of threads to allocate to the blast call (default is 1)

* **-op, --only_prepare**: creates the blast databases, but only prints blast commands to run on STDOUT). 

## 1.4.1 HOGs2clus.py

This script converts the phylogenetic hierarchical orthogroups (HOGs) created by [OrthoFinder](https://github.com/davidemms/OrthoFinder) >= 2.4 into clus format.

There are two positional arguments, results are saved in the specified output file:

usage: `HOGs2clus.py input output`

* **input**: name of the HOG file to convert, tab-separated values (e.g N0.tsv, located in the "Phylogenetic Hierarchical Orthogroups" directory in Orthofinde results)

* **output**: name of the output clus file

## 1.4.1 cleanup_clus.py

This script cleans up a clus file. Since HOGs do not include singletons (genes unassigned to orthogroups), this scripts can add singleton orthogroups. This script can also be used to only keep a clus file with sequences of interest.

There are two positional arguments, the  results are printed to STDOUT.

usage: `cleanup_clus.py acc_list clus_file > output.clus`

+ **acc_list**: newline-delimited list of accessions that should be in the output clus file. Accessions in the input and acc list will be kept, others will be deleted. Accessions not in the input clus file will be added as singletons

+ **clus_file**: name of the input clus file to clean



## 1.5 prepare_homo_MCL_length_file.py

For coarse-grained orthology, one might wanna use [homomcl](https://github.com/willpett/homomcl). That requires a specific length file. This script creates such a file from a fasta file.

usage: `prepare_homomcl_length_file.py proteome > lengthfile`
