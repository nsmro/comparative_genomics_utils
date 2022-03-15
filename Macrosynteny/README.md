These are some scripts for doing oxford dotplots

# 1. List of the tools

The scripts work with **msynt files** (see description of the format in [here](#msynt))

## 1.1 make_chromosome_relabel.py

Makes a tidy table mapping chromosomes/scaffold names to length-based names (longer  chromosome gets names prefix01, second longest prefix 02, etc.).

The following arguments are required:

* **genome (positional)**: a genome file, fasta format

* **-p, --prefix**: prefix to add to chromosome numbers (e.g. for *Euprymna scolopes*, coule be EUPSC_, Eupsc, or Esc)

* **-n, --number chromosomes**: how many of the largest chromosomes should be kept. any number higher than the number of scaffolds will lead to all the scaffolds to be renamed

This prints to STDOUT three fields: new name, old name and length of the chromosome/scaffold (in basepairs)

## 1.2 relabel_msynt.py

This script will glob all files with the *msynt* extension located in the current working directory, and will write a renamed file for each msynt file it globs (adding the "relab" infix before the "msynt suffix". For example for a `EUPSC.msynt` input file, it will write a `EUPSC.relab.msynt` file).

Usage : `relabel_msynt.py mylabels.labels prefixlength`

The two positional arguments that are required are:

* **maylabels.labels**: a file generated using `make_chromosome_relabel.py` 

* **prefixlength**: the number of characters in the specified prefix ("EUPSC_" is of length 6, "Eupsc" is of length 5, and "Esc" of length 3)

This is done to relabel synt files created beforehand with original scaffold. For example, if you wish your the chromosome on your dotplots to have a different name than the ones in the genome.

## 1.3 macrosynteny_homology_dotplot.R

This script requires the following R packages: gtools, slanter, RColorBrewer, tidyverse and argparse.

This script requires at least two msynt files. And for every pair of msyntfiles, it'll plot an oxford dotplot of the chromosomes. Orthology is stated by in msynt file. Only relative positions will be plotted. Each dot of the oxford dotplot corresponds to an inter-species orthology relationship (cartesian product of the orthogroup members of the two species)

The parameters are:

* **msynt (positinal)**: space-delimited list of msyntfiles. At least two files. Filter which chromosomes you wwant to keep before runninng this script. The species prefix should be the first field of the name (dot-delimited, e.g. the prefix "EUIPSC" will be isolated from the filename "EUPSC.syntenyrun.relabel.mysynt")

* **--output_folder**: name of the folder where the ouput files should be returned. Default: current working directory.

* **--max_para**: for a given animal maximum number of paralogs an OG is allowed to have

* **--min_genes**: minimum number of disting orthologs a chromosome/scaffold should bear to be retained

* **--reorder**: if this flag is not selected, the chromosomes will be ordered alphabetically along the x and y axis (for a clean order, the first chromosome of an animal with the prefix Esc should rather be called Esc01 instead of Esc1). If this flag is reordered, the rows and columns will be clusterd using the Ward's minimum variance method, and the `slanter` package will be used to align the highest number of shared orthologs onto the diagonal.

For every pair of msyntfiles  "SPECIES1.msynt" and "SPECIES2.msynt", assuming the reorder flag has not been provided and `--maxpara` is of 1 three pdf files will be generated:

- `SPECIESX-SPECIESY.maxpara_1.dotplot.pdf`: the oxford dotplot

- `SPECIESX-SPECIESY.maxpara_1.bubble_plot.fisher_enrichments.pdf`: a significance grid, where circle size corresponds to the number of orthologues a pair of chromosome share, and the color the pvalue

- `SPECIESX-SPECIESY.maxpara_1.dotplot.pdf`: the heatmap with the number of shared genes.



<a id="msynt"></a>

# 2. description of the *msynt* format

The msynt format is used only for the aforementionned scripts to run. The msynt format comprises five tab-separated fields.

The five fields are as follows:

1. **scaffold**: name of the scaffold/schromosome on which the gene is located

2. **accession**: accession of the gene

3. **position**: start position of the gene on the scaffold (the start position is used to define the indices of the genes)

4. **orthogroup**: id of the orthogroup/orthogroup name. For the macrosynteny R script, the id should be consistent across the msynt files compared

5. **nb_para**: number of paralogues of the orthogroup belong to the species


