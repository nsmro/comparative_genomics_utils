#!/usr/bin/env Rscript

library(gtools)
library(slanter)
library(RColorBrewer)
library(tidyverse)
library(argparse)



parser <- ArgumentParser(description = 'Provided a list of msynt files, plots a bubble plot showing homologies, as well as dotplot for each msynt file pair.')
parser$add_argument('msynt', nargs = '+', type = 'character', help = 'space-separated list of msynt files. filter chromosomes you want to keep before running this script.')
parser$add_argument('--output_folder',  type = 'character', help = 'name of the output where all the files will be created. If not specified, everything will be saved in current working dir', default = getwd())
parser$add_argument('--max_para',  type = 'integer', help = 'max number of paralogs retained genes are allowed to have', default = 4)
parser$add_argument('--min_genes',  type = 'integer', help = 'minimum of distinct orthogroups a scaffold should carry to be retained', default = 2)
parser$add_argument('--reorder',  help = 'use this flag if you want to reorder chromosomes based on shared gene clustering', action = "store_true", default = FALSE)
args <- parser$parse_args()



#' load msynt files 
#'
#' @param msyntlist a vector of msynt filenames
#'
#' @return list of msynt tibbles, each can be accessed by prefix name
load_msyntfiles <- function(msyntlist){
    cols_msynt <- c('scaffold', 'accession', 'position', 'orthogroup', 'nb_para')
    output <- list()
    for (msyntfile in msyntlist){
        prefix <-  strsplit(x = basename(msyntfile),
                          split = ".",
                          fixed = T)[[1]][1]
        tmp_list <- readr::read_delim(file = msyntfile,
                                        col_names = cols_msynt,
                                        delim = "\t") %>% list(.)# wrap the tibble inside a list before adding it to output list
        names(tmp_list) <- prefix
        output <- c(output, tmp_list)
    }
    return(output)
}

#' This is a function for applying mixedorder on a tibble using dplyr::arrange
mixedrank <- function(x) order(gtools::mixedorder(x))

#' filters a msynt tibble
#'
#' @param msyntlist a list of loaded *.msynt files
#' @param prefix a prefix filename used to pick up the raw syntfile from msyntlist
#' @param max_para max number of paralogs retained genes are allowed to have
#' @param min_genes minimum of distinct orthogroups a scaffold should carry to be retained
#'
#' @return tibble with 7 columns (scaffold,accession,position,orthogroup,nb_para,nb_gene_scf,nb_unique_gene_fam_scf) sorted by position in bp
process_msynt <- function(msyntlist,
                          prefix,
                          max_para = args$max_para,
                          min_genes = args$min_genes){
    if (identical(args$microsynt_labelling, F) == F){ #force keep microsyntenic genes
        species_lab <- ls_microsynt_labelling[[prefix]] %>% dplyr::rename(synt_orthogroup = orthogroup)
        tmp_output <- msyntlist[[prefix]]
        tmp_output <- tmp_output %>%
                        dplyr::left_join(species_lab, by = 'accession') %>%
                        dplyr::filter(synt_orthogroup != NA | nb_para <= max_para)
        tmp_output$orthogroup <- dplyr::coalesce(tmp_output$synt_orthogroup, tmp_output$orthogroup) #correspondance between syntenic OGs only
 
    } else{
        tmp_output <- msyntlist[[prefix]]  %>%
                    dplyr::filter(nb_para <= max_para)
    }
        output <- tmp_output %>% dplyr::group_by(scaffold) %>%
                    dplyr::mutate(nb_gene_scf = n(), # genes by scaffold
                                  nb_unique_gene_fam_scf = length(unique(orthogroup))) %>% # unique OGs by scaffold
                    dplyr::ungroup() %>%
                    dplyr::filter(nb_unique_gene_fam_scf >= min_genes) %>%
                    dplyr::arrange(position) %>%
                    dplyr::arrange(mixedrank(scaffold)) %>%
                    dplyr::select(c('scaffold','accession','position','orthogroup','nb_para','nb_gene_scf','nb_unique_gene_fam_scf'))
    return(output)
}


#' filters msyntfile based on a user-provided list of chromosomes, sort the chroms by name and position
#' adds indices based on relative position (e.g. genes on chrom 1: 1-10, on chrom2: 11-15 etc.)
#' calculate location of chromosome breaks
#' @param tbl a tibble, output of process_msynt
#' @param chromlist vector of chromosomes to filter the table
#'
#' @return filtered tibble (based on.c) with 7 columns (scaffold,accession,position,orthogroup,nb_para,nb_gene_scf,nb_unique_gene_fam_scf) sorted by position in bp
#' @return chromosome breaks indices
#' @return chromosome names
add_indices_get_bks <- function(tbl){
    ouput_tbl <- tbl %>%
             dplyr::mutate(id = row_number())
    tmp_tbl <- dplyr::count(ouput_tbl, scaffold) %>%
             dplyr::arrange(mixedrank(scaffold))
    chrom_bks <- cumsum(c(1, tmp_tbl$n))
    chrom_names <- unique(tmp_tbl$scaffold)
    output <- list(ouput_tbl, chrom_bks, chrom_names)
    names(output) <- c('sorted_msynt', 'chromosome_breaks', 'chromosome_names')
    return(output)
}


#' makes matrices of shared genes for heatmap
#'
#' @param msyntfile_a a tibble, output of process_msynt (columns)
#' @param msyntfile_b a tibble, output of process_msynt
#'
#' @return list of two matrices, one with normalized values, one with absolute number of shared genes
matrices_shared_genes <- function(msynt_a, msynt_b){
    chromlist_a <- unique(msynt_a$scaffold)
    chromlist_b <- unique(msynt_b$scaffold)
    output_matrix <- rep(list(matrix(ncol = length(chromlist_a),
                                     nrow = length(chromlist_b),
                                     dimnames = list(chromlist_b, chromlist_a))),
                         2)
    names(output_matrix) <- c("norm", "abs")
    for (chrom_i in chromlist_a){
        for(chrom_j in chromlist_b){
            fam_in_i <- unique(dplyr::filter(msynt_a, scaffold == chrom_i)$orthogroup) #nb OGs of chrom_i of sp_x
            fam_in_j <- unique(dplyr::filter(msynt_b, scaffold == chrom_j)$orthogroup) #nb OGs of chrom_j of sp_y
            max_shared <- min(length(fam_in_i), length(fam_in_j)) #max nb of possible shared families between i and j (for norm)
            shared_ogs <- sum(fam_in_i %in% fam_in_j)  #calculate shared number of families
            output_matrix$norm[chrom_j, chrom_i] <- shared_ogs / max_shared
            output_matrix$abs[chrom_j, chrom_i] <- shared_ogs
        }
    }
return(output_matrix)
}

#This will be a list of tibble, one item appended to it each iteration
gene_enrichment_table <- list()

msyntfile_ls <- load_msyntfiles(args$msynt)

for (file_pair in combn(msyntfile_ls, 2, simplify = F)){
    species_X <- names(file_pair[1])
    species_Y <- names(file_pair[2])

    msynt_x <- process_msynt(msyntfile_ls, species_X)
    msynt_y <- process_msynt(msyntfile_ls, species_Y)

    #all this part up to bbplot is for making heatmaps and determining chrom homology
    matrix_x_y <- matrices_shared_genes(msynt_x, msynt_y)
    if (args$reorder == T){reorder <- '.reordered'} else{reorder <- ''}

    output_prefix <- paste0(args$output_folder, '/', species_X,'-',species_Y, '.maxpara_', args$max_para)
    output_heatmap_current_pair <- paste0(output_prefix, '.heatmap.shared.genes', reorder, '.pdf')
    output_bubble_current_pair <- paste0(output_prefix, '.bubble_plot.fisher_enrichments', reorder, '.pdf')
    dotplot_output_filename <- paste0(output_prefix, '.dotplot', reorder, '.pdf')
    
    #Now for the heatmap, clustering is done on both rows and columns
    heatmap_norm <- slanter::sheatmap(matrix_x_y$norm,
             color = colorRampPalette(RColorBrewer::brewer.pal(n = 7, name = 'YlOrRd'))(100),
             cluster_rows = T,
             cluster_cols = T,
             clustering_distance_rows = 'euclidean',
             clustering_distance_cols = 'euclidean',
             clustering_method = 'ward.D2',
             filename = output_heatmap_current_pair)
    custom_order <- slanter::slanted_orders(matrix_x_y$norm, order_rows = TRUE, order_cols = TRUE)
    #However, the order from the clustering is only used if the "reorder" flag is used in arguments.
    if (args$reorder == T){
        row_labels <- rownames(matrix_x_y$abs[custom_order$rows, ])
        col_labels <- colnames(matrix_x_y$abs[,custom_order$cols])
    } else{
        row_labels <- gtools::mixedsort(rownames(matrix_x_y$norm))
        col_labels <- gtools::mixedsort(colnames(matrix_x_y$norm))
    }

    #sanity check
    m_abs_sorted <- matrix_x_y$abs[row_labels, col_labels]
    m_norm_sorted <- matrix_x_y$norm[row_labels, col_labels]
    
    # Fisher enrichment test, stored in a matrix
    fisher_pvalues <- matrix(ncol = length(colnames(m_abs_sorted)),
                             nrow = length(rownames(m_abs_sorted)),
                             dimnames = dimnames(m_abs_sorted))
    
    for (chr_y in rownames(m_abs_sorted)){
        for (chr_x in colnames(m_abs_sorted)){
            not_x <- colnames(m_abs_sorted) != chr_x
            not_y <- rownames(m_abs_sorted) != chr_y
            contingency_x_y <- matrix(c(m_abs_sorted[chr_y,chr_x],
                                        sum(m_abs_sorted[not_y,chr_x]),
                                        sum(m_abs_sorted[chr_y,not_x]),
                                        sum(m_abs_sorted[not_y,not_x])),
                                      nrow=2,ncol=2)
            fisher_pvalues[chr_y, chr_x] <- fisher.test(contingency_x_y, alternative = 'greater')$p.value
        }
    }
    
    #correct for multiple tests with Benjamini-Hochberg
    fisher_pvalues_fdr <- matrix(p.adjust(as.vector(fisher_pvalues), method = "BH"),
                                 ncol = length(colnames(m_abs_sorted)),
                                 nrow = length(rownames(m_abs_sorted)),
                                 dimnames = dimnames(m_abs_sorted))
    

    #make a dataframe for bubble plot data.frame parser already converts contingency tables (obtained via as.table) to long dataframe
    df_abs_sorted <- tibble::tibble(data.frame(as.table(m_abs_sorted)),
                                    .name_repair = ~ c(species_Y, species_X, 'shared_genes'))
    df_fisher_pvalues_fdr <- tibble::tibble(data.frame(as.table(fisher_pvalues_fdr)),
                                    .name_repair = ~ c(species_Y, species_X, 'pval'))
    bbplot_table <- dplyr::full_join(df_abs_sorted, df_fisher_pvalues_fdr, by = c(species_X, species_Y))

    bbplot_table[species_X] <- factor(dplyr::pull(bbplot_table, species_X), #contrary to `$`syntax, dplyr::pull() allows to use string stored in a variable
                                      levels = col_labels)
    bbplot_table[species_Y] <- factor(dplyr::pull(bbplot_table, species_Y),
                                      levels = row_labels)

    #This is just to make the ggplot call less cluttered
    bbplot_table$pval_factor <- cut(bbplot_table$pval,
                                    breaks = c(0, 0.1, 0.5, 1),
                                    labels = c("<0.01","<0.05",">0.05"),
                                    include.lowest = T)
    
    bks_shared_genes_size <- c(0,50,100,150) 
    bks_shared_genes_label <- c(">=1",">=50",">=100",">=150")
    g <- ggplot2::ggplot(bbplot_table, aes_string(x = species_X, y = species_Y, size = 'shared_genes', color = 'pval_factor')) +
          ggplot2::geom_point(alpha = 0.5) +
          ggplot2::scale_size(name = "number of shared genes",
                              breaks = bks_shared_genes_size,
                              labels = bks_shared_genes_label,
                              guide = "legend")+
          ggplot2::scale_color_manual(breaks = c("<0.01","<0.05",">0.05"),
                                      values = c("red","orange", "grey"),
                                      name = "p-values",
                                      drop = F) + #even if a pvalue threshold isn't populated, have it displayed in the legend
          ggplot2::theme_minimal() +  
          theme(axis.title.x = element_text(size = 8,
                            hjust = 0.95,
                            vjust=0.5),
                axis.title.y = element_text(size = 8,
                                            hjust = 0.95,
                                            vjust=0.5),

                axis.text.x = element_text(angle = 90,
                                           size = 6,
                                           hjust = 0.95,
                                           vjust=0.5),
                axis.text.y = element_text(size = 6,
                                           hjust = 0.95,
                                           vjust=0.5))
    
        ggsave(plot = g,
           file = output_bubble_current_pair,
           unit = 'cm',
           width = 15,
           height = 12)


    #saving chromosome pairs info into the list initiated before the for loop
    pvalues_tbl <- data.frame(fisher_pvalues_fdr) %>%
                        tibble::rownames_to_column() %>%
                        tidyr::gather(key = 'key', value = 'value', -rowname)
    names(pvalues_tbl) <- c('Y', 'X', 'pvalue')
    norm_counts_tbl <- data.frame(m_norm_sorted) %>%
                        tibble::rownames_to_column() %>%
                        tidyr::gather(key = 'key', value = 'value', -rowname)
    names(norm_counts_tbl) <- c('Y', 'X', 'norm_gene_count')
    raw_counts_tbl <- data.frame(m_abs_sorted) %>%
                        tibble::rownames_to_column() %>%
                        tidyr::gather(key = 'key', value = 'value', -rowname)
    names(raw_counts_tbl) <- c('Y', 'X', 'gene_count')
    pair_info <- pvalues_tbl %>%
                        dplyr::inner_join(norm_counts_tbl, by = c('Y', 'X'))%>%
                        dplyr::inner_join(raw_counts_tbl, by = c('Y', 'X'))
    pair_info$X <- paste(species_X, pair_info$X, sep = '_')
    pair_info$Y <- paste(species_Y, pair_info$Y, sep = '_')
    gene_enrichment_table <- c(gene_enrichment_table, list(pair_info))
    

    #from now on, this is for making macrosynteny dotplots
    #we reorder rows based on order defined previously
    #add indices and get chromosome breaks for each chrom
    
    fam_in_x <- unique(msynt_x$orthogroup) #OGs in sp_x
    fam_in_y <- unique(msynt_y$orthogroup) #OGs in sp_y
    shared_ogs <- intersect(fam_in_y, fam_in_x)  #calculate shared families


    msynt_x$scaffold <- factor(msynt_x$scaffold,  levels = col_labels)
    msynt_x <- msynt_x %>% dplyr::filter(orthogroup %in% shared_ogs) %>%
                    dplyr::arrange(key = scaffold)
    tmp_ls_x <-  add_indices_get_bks(msynt_x) # add indices on the filtered msynt file
    msynt_x <- tmp_ls_x$sorted_msynt
    chrom_bks_x <- tmp_ls_x$chromosome_breaks
    chrom_names_on_x <- tmp_ls_x$chromosome_names

    msynt_y$scaffold <- factor(msynt_y$scaffold,  levels = row_labels)
    msynt_y <- msynt_y %>% dplyr::filter(orthogroup %in% shared_ogs) %>%
                    dplyr::arrange(key = scaffold)
    tmp_ls_y <-  add_indices_get_bks(msynt_y)

    msynt_y <- tmp_ls_y$sorted_msynt
    chrom_bks_y <- tmp_ls_y$chromosome_breaks
    chrom_names_on_y <- tmp_ls_y$chromosome_names


    #rename cols to keep control on their names
    rename_list_x <- c('scaffold', 'accession', 'id')
    names(rename_list_x ) <- c(paste(species_X, 'scaffold', sep = '_'), paste(species_X, 'accession', sep = '_'),  paste(species_X, 'id', sep = '_'))
    msynt_x_renamed <- msynt_x[c('scaffold', 'accession', 'orthogroup', 'id')] %>% #select only cols that are relevant
                           dplyr::rename(all_of(rename_list_x)) #rename columns using named vector where new/old names are vector's names/values
    
    rename_list_y <- c('scaffold', 'accession', 'id')
    names(rename_list_y ) <- c(paste(species_Y, 'scaffold', sep = '_'), paste(species_Y, 'accession', sep = '_'),  paste(species_Y, 'id', sep = '_'))
    msynt_y_renamed <- msynt_y[c('scaffold', 'accession', 'orthogroup', 'id')] %>% #select only cols that are relevant
                           dplyr::rename(all_of(rename_list_y)) #rename columns using named vector where new/old names are vector's names/values
    
    
    dot_plot_table <- msynt_x_renamed %>%
                        dplyr::inner_join(msynt_y_renamed, by = 'orthogroup') #Ortho found in both species, genes of one are repeated if multiple matches
    dot_plot_table$alpha <- 0.8


    xmax <- max(msynt_x$id) #The max indices and the ones of the breaks are the unique gene families.
    ymax <- max(msynt_y$id) #If one indice is absent of one species, it's dropped, and if extremity is dropped, this value is messed up
    
    
    #Position in the middle of consecutive breaks
    chrom_position_on_x <- (chrom_bks_x[- c(1)] + chrom_bks_x[- c(length(chrom_bks_x))]) / 2
    chrom_position_on_y <- (chrom_bks_y[- c(1)] + chrom_bks_y[- c(length(chrom_bks_y))]) / 2
    
    #if no labelling specified, is_null(df_labels) returns TRUE

    myggproto <- ggplot2::ggplot(dot_plot_table, aes_string(x = paste(species_X, 'id', sep = '_'),
                                                            y = paste(species_Y, 'id', sep = '_'),
                                                            alpha = 'alpha'))
    
    g <- myggproto +
            ggplot2::geom_point(size = 0.8, stroke = 0) + #stroke = 0 avoids ringing effect when using alpha
            ggplot2::theme_linedraw() +
            ggplot2::scale_alpha(guide = 'none') +
            ggplot2::scale_y_continuous(breaks = chrom_bks_y,
                                        sec.axis = dup_axis(labels = NULL)) +
            ggplot2::scale_x_continuous(breaks = chrom_bks_x,
                                        sec.axis = dup_axis(labels = NULL)) +
            ggplot2::theme(plot.margin = unit(c(3,3,1,1), "cm"),# top,right,left,bottom
                           axis.title.x = element_text(size = 8,
                                                       hjust = 0.95,
                                                       vjust=0.5),
                           axis.title.y = element_text(size = 8,
                                                       hjust = 0.95,
                                                       vjust=0.5),
                           axis.title.y.right = element_blank(),
                           axis.title.x.top = element_blank(),
                           axis.text.y = element_text(size = 6),
                           axis.text.x = element_text(size = 6, angle = 90),
                           panel.grid.major = element_line(colour = "black",
                                                           linetype = "dashed",
                                                           size = 0.1), #lines for chromosome breaks
                           panel.grid.minor = element_blank(),
                           axis.ticks.length.x.top = unit(1.5, "cm"),
                           axis.ticks.length.y.right = unit(1.5, "cm"),
                           legend.position = 'bottom') +
            ggplot2::labs(y = species_Y, x = species_X) +
            ggplot2::annotate("text",
                              x = xmax + 100,
                              y = chrom_position_on_y,
                              hjust = 0,
                              label = chrom_names_on_y,
                              size = 1.7) + # 
            ggplot2::annotate("text",
                              x = chrom_position_on_x,
                              y = ymax + 100 ,
                              hjust = 0,
                              angle = 90,
                              label = chrom_names_on_x,
                              size = 1.7) +
            ggplot2::coord_cartesian(xlim = c(1, xmax + 2),
                                     ylim = c(1, ymax + 2),
                                     expand = FALSE,
                                     clip = "off") # chrom labels are outside the plotting area, we don't want to clip them
    
    ggsave(plot = g,
           file = dotplot_output_filename,
           unit = 'cm',
           width = 15,
           height = 18)

}

#merge the list of tibbles into a single one, and write the the output file
gene_enrichment_table <- dplyr::bind_rows(gene_enrichment_table)
readr::write_delim(x = gene_enrichment_table,
                   path = paste0(args$output_folder, '/gene_enrichment_table.tsv'),
                   col_names = FALSE,
                   delim = '\t')
