import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from plotting_utils import reading_file

from plotting_utils import mm2inch, rcParams
for k, v in rcParams.items():
    plt.rcParams[k] = v

# Loading HGNC translation
hgnc_df = pd.read_csv(snakemake.input.HGNC, sep='\t')
hgnc_df['ncbi_id'] = "GeneID:" + hgnc_df['ncbi_id'].astype('str').str.split('.').str[0]

# Creating translation dictionary
ref2hgnc = hgnc_df.set_index('ncbi_id')['HGNC_ID'].to_dict()
ens2hgnc = hgnc_df.set_index('ensembl_id')['HGNC_ID'].to_dict()
gene2hgnc = {**ref2hgnc, **ens2hgnc}

# Loading genes
up, down = reading_file(snakemake.input.gene_list)
up = [ens2hgnc[gene] for gene in up]
down = [ens2hgnc[gene] for gene in down]

# Loading data
overlaps = dict()
for overlap in snakemake.input.overlaps:
    df = pd.read_csv(overlap, sep='\t', names=['annotation', 'gene', 'score'])
    annotation, gene, score = df.values.tolist()[0]

    if annotation not in overlaps.keys():
        overlaps[annotation] = dict()
    overlaps[annotation][gene2hgnc[gene]] = score

data = pd.DataFrame.from_records(data=overlaps)
data.fillna(value=0, inplace=True)

# Plotting
fig, axes = plt.subplots(ncols=2, nrows=1, figsize=mm2inch((70,50)))
data *= 100
for ax, ori in zip([0, 1], [down, up]):
    sns.violinplot(x='variable', y='value', data=data[data.index.isin(ori)].melt(), cut=0, bw=0.25, ax=axes[ax])
    axes[ax].set_ylim([0, 100])
axes[1].get_yaxis().set_visible(False)
axes[0].set_ylabel('% of exon overlapped')
axes[0].set_title('EM'+snakemake.wildcards.comp+' down')
axes[1].set_title('EM'+snakemake.wildcards.comp+' up')

fig.subplots_adjust(wspace=0, right=0.99, bottom=0.10, top=0.90, left=0.10)
plt.savefig(snakemake.output.plot)
