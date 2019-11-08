import itertools
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

from get_gene_gtf import get_gene

df_gene, gene_start, gene_end, gene_strand, gene_seq = get_gene(
    snakemake.input.gtf, snakemake.wildcards.gene
)

# Extracting information about one transcript
transcript = sorted(df_gene['transcript_name'].dropna().unique().tolist())[0]
df_transcript = df_gene[df_gene['transcript_name'] == transcript]
df_transcript = df_transcript[df_transcript['feature'] == 'exon']

exons = df_transcript[['start', 'end']].values
exon = exons[4]

# Loading gene coverage
data = pd.read_csv(snakemake.input.gene_coverage, sep='\t', index_col=0, header=[0, 1, 2, 3, 4])

# Removing annotation since only STAR uses a GTF
data = data.groupby(level=['aligner', 'dataset', 'tissue', 'trimmer'], axis=1).mean()
data = data.iloc[int(exon[0])-gene_start-2:int(exon[1])-gene_start+1]

# Plotting, for a tissue at a time?
fig, axes = plt.subplots(nrows=2, ncols=2, sharey=True, sharex=True)
tissues = sorted(list(set(data.columns.get_level_values('tissue'))))
perms = [(0, 0), (0, 1), (1, 0), (1, 1)]

for tissue, indx in zip(tissues, perms):

    col_tissue = data.columns[data.columns.get_level_values('tissue') == tissue]
    data_tissue = data[col_tissue]

    # For aligners
    df_mean = data_tissue.groupby(level='aligner', axis=1).mean()
    df_std = data_tissue.groupby(level='aligner', axis=1).std()

    for x in range(len(df_mean.columns)):
        axes[indx].plot(df_mean.index, df_mean[df_mean.columns[x]], linewidth=4)

        axes[indx].fill_between(
            x = df_mean.index,
            y1 = df_mean[df_mean.columns[x]] + df_std[df_std.columns[x]],
            y2 = df_mean[df_mean.columns[x]] - df_std[df_std.columns[x]],
            alpha=0.3
        )

        axes[indx].set_title(tissue)

plt.savefig(snakemake.output.plot)
