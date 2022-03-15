#!/usr/bin/env python3

import argparse
import scipy
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

parser = argparse.ArgumentParser(description = """Makes a heatmap out of the provided trancript table (rows = transcripts, columns = stages)""")
parser.add_argument('-b', '--bar_label', help = 'Caption of the legend bar.', default = 'zscore')
parser.add_argument('input', help = 'input table. Columns are stages, rows are transcripts')
parser.add_argument('-s', '--title', help = 'title of the plot. Default is args.input.')
parser.add_argument('-o', '--output', help = 'name of the output file, default name is as follows: if input is "inputname.tmps.tsv" output is "inputname_heatmap.pdf')
parser.add_argument('-t', '--threshold', help = 'Mask tpms under certain threshold.', type = float, default = 1)
parser.add_argument('-c', '--clustermap', help = 'If the flag is used, rows will be clustered together', action = 'store_true')
parser.add_argument('-r', '--rectangles', help = 'If the flag is used, each field of the heatmap a rectangle (default: square). In the case of clustermap, no squares', action = 'store_false')
parser.add_argument('--width', help = 'figure width, in inches, default is 8', default = 8, type = int)
parser.add_argument('--height', help = 'figure height, in inches, default is 6', default = 6, type = int)
args = parser.parse_args()


if args.title is None:
    args.title = args.input
if args.output is None:
    output_prefix = args.input.split('.')[0]
    args.output = f'{output_prefix}_{args.threshold}masked_heatmap.pdf'

print(f'Heatmap will be saved under {args.output}')

df = pd.read_csv(args.input, index_col=0, sep = '\t')
df = df.dropna(axis = 1)

df_mask = df.apply(lambda x: x < args.threshold)

z_score_df = df.apply(scipy.stats.mstats.zscore, axis = 1, result_type = 'expand')

#if the zscore is NaN, means it isn't expressed at all
is_NaN = z_score_df.isnull()
row_has_NaN = is_NaN.any(axis=1)
ls_na = z_score_df[row_has_NaN].index
string_ls_na = ', '.join(ls_na)
print(f'Dropped accessions (not expressed): {string_ls_na}')

z_score_df = z_score_df.dropna(axis = 0)
df_mask = df_mask.drop(index = ls_na) #correct the mask to delete values not in zscoredf

# apply zscore normalization to every row (axis = 1)
#result_type expand makes an output into columns. If not, we end up with 2 columns, one with names, one with the array of results
z_score_df.columns = list(df.columns.values)

if args.clustermap is False:
    fig,ax = plt.subplots(figsize = (args.width, args.height))
    ax = sns.heatmap(data = z_score_df,
                    ax = ax,
                    linewidths = .15,
                    square = args.rectangles,
                    center = 0,
                    cmap = 'RdYlBu_r',
                    mask = df_mask,
                    vmin = -4,
                    vmax = 4, # maximum values on the colorbar
                    cbar_kws={'label': args.bar_label, 'shrink': .3, 'ticks':list(range(-4,5,2))})# dictionary for colorbar. shrink it so that it doesn't look too big next to the heatmap
    ax.set_facecolor('#CCCCCC')
    bottom, top = ax.get_ylim()
    #ax.set_ylim(bottom + 0.5, top - 0.5) #Matplotlib broke seaborn heatmaps... So wating for it to be fixed, here is a little bit of code to fix...
    ax.set_title(f'{args.title}', fontsize = 10)
    ax.set_xticklabels(ax.get_xmajorticklabels(), fontsize = 5)
    ax.xaxis.set_ticks_position('none') 
    ax.set_yticklabels(ax.get_ymajorticklabels(), fontsize = 5)
    ax.yaxis.set_ticks_position('none')
    plt.tight_layout()
    plt.savefig(args.output)
else:
    g = sns.clustermap(data = z_score_df,
                    linewidths = .15,
                    square = False,
                    center = 0,
                    cmap = 'RdYlBu_r',
                    vmin = -4,
                    vmax = 4, # maximum values on the colorbar
                    mask = df_mask,
                    col_cluster = False,
                    yticklabels = True,
                    figsize = (args.width, args.height),
                    cbar_kws={'label': args.bar_label, 'shrink': .3, 'ticks':list(range(-4,5,2))})# dictionary for colorbar. shrink it so that it doesn't look too big next to the heatmap
    g.ax_heatmap.set_facecolor('#CCCCCC')
    g.ax_heatmap.set_title(f'{args.title}', fontsize = 10)
    g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xmajorticklabels(), fontsize = 5)
    g.ax_heatmap.xaxis.set_ticks_position('none') 
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_ymajorticklabels(), fontsize = 5)
    g.ax_heatmap.yaxis.set_ticks_position('none') 
    g.savefig(args.output)
