import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd

from snakemake import shell

from plotting_utils import mm2inch, rcParams
for k, v in rcParams.items():
    plt.rcParams[k] = v


def getStrand(key):
    for annotation in annotations:
        if key in data[annotation].keys():
            strand = data[annotation][key]['strand'].values[0]

    if (key == HGNC2ens[HGNC]) or (key == HGNC2ref[HGNC]):
        return 1
    elif strand == '+':
        return 2
    elif strand == '-':
        return 0


# Loading parameters
HGNC = snakemake.wildcards.HGNCgene
annotations = ['ensembl92', 'ensembl98', 'refseq']
columns = ['chr', 'start', 'end', 'gene', 'score', 'strand']

# Loading HGNC
hgnc_df = pd.read_csv(snakemake.input.HGNC, sep='\t')
hgnc_df['ncbi_id'] = "GeneID:" + hgnc_df['ncbi_id'].astype('str').str.split('.').str[0]
HGNC2ens = hgnc_df.set_index('HGNC_ID')['ensembl_id'].to_dict()
HGNC2ref = hgnc_df.set_index('HGNC_ID')['ncbi_id'].to_dict()
ens2HGNC = hgnc_df.set_index('ensembl_id')['HGNC_ID'].to_dict()
ref2HGNC = hgnc_df.set_index('ncbi_id')['HGNC_ID'].to_dict()
gene2HGNC = {**ens2HGNC, **ref2HGNC}

# Figure params
line_weight = 0.025

# Loading beds
data = dict()
max, min = 0, 1e12
for annotation in annotations:
    data[annotation] = dict()
for bed in snakemake.input.beds:
    df = pd.read_csv(bed, sep='\t', names=columns)
    df['HGNC'] = df['gene'].map(gene2HGNC)
    gene = df['gene'].values[0]
    for annotation in annotations:
        if annotation in bed:
            data[annotation][gene] = df
            _max = np.max(np.array(df[['start', 'end']]))
            _min = np.min(np.array(df[['start', 'end']]))
            if _max > max:
                max = _max
            if _min < min:
                min = _min


# Plotting?
ratios = [len(data[annotation]) for annotation in annotations]
fig, axes = plt.subplots(ncols=1, nrows=3, gridspec_kw={'height_ratios':ratios}, figsize=mm2inch((178, 50)))

for ax, annotation in enumerate(annotations):
    n_gene = len(data[annotation])

    axes[ax].set_ylim([0, n_gene])
    axes[ax].set_xlim([min-1000, max+1000])

    for idx, gene in enumerate(sorted(data[annotation], key=getStrand)):
        color = 'k'
        if (gene == HGNC2ens[HGNC]) or (gene == HGNC2ref[HGNC]):
            color = 'g'

        gene_max, gene_min = 0, 1e12
        for _, row in data[annotation][gene].iterrows():
            if row.start < gene_min:
                gene_min = row.start
            if row.end > gene_max:
                gene_max = row.end

            axes[ax].add_patch(Rectangle(
                (row.start, idx+0.25), row.end-row.start, 0.5, color=color
            ))

        axes[ax].add_patch(Rectangle(
            (gene_min, idx+0.25+(0.5/2)-(line_weight/2)),
            gene_max-gene_min, line_weight, color=color
        ))

        for x in np.arange(gene_min, gene_max, (max-min)/80)[1:]:
            if row.strand == '+':
                x = [x, x+(max-min)/80/5, x]
            else:
                x = [x, x-(max-min)/80/5, x]

            y = [idx+0.5-line_weight*3, idx+0.5, idx+0.5+line_weight*3]
            axes[ax].add_line(Line2D(
                x, y, lw=1., color=color
            ))

plt.tight_layout()
plt.subplots_adjust(hspace=0)
plt.savefig(snakemake.output.plot)
