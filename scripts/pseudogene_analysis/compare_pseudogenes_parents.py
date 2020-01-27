import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

from plotting_utils import mm2inch, rcParams
for k, v in rcParams.items():
    plt.rcParams[k] = v


def reading_file(path):
    genes = list()
    with open(path, 'r') as f:
        for line in f.readlines():
            line = line.strip()
            genes.append(line)
    return genes

# Reading files
# Using data from http://pseudogene.org/psicube/index.html
df_gene_biotype = pd.read_csv(snakemake.input.gene_biotype, sep='\t')
df_parent = pd.read_csv(snakemake.input.parent, sep='\t')
df_biotype = pd.read_csv(
    snakemake.input.biotype, sep='\t',
    names=[
        'seqname', 'source', 'feature', 'start',
        'end', 'score', 'strand', 'frame', 'attributes'
    ]
)

# Reading data genes
data = pd.read_csv(
    snakemake.input.data, sep='\t',
    index_col=[0, 1, 2, 3], header=[0, 1, 2, 3, 4, 5],
)
data.dropna(axis=0, inplace=True)

# Reading gene list
lgenes = reading_file(snakemake.params.gene_list)
lgenes2 = reading_file(snakemake.params.gene_list2)

# Parsing df_biotype
extracting = {
    "gene_id": "gene_id \"(.*?)\"",
    "biotype": "transcript_type \"(.*?)\"",
    "transcript_id": "transcript_id \"(.*?)\""
}
for k, re in extracting.items():
    df_biotype[k] = df_biotype['attributes'].str.extract(re)

# Keeping only processed_pseudogenes and transcribed_processed_pseudogenes
accepted = ['processed_pseudogene', 'transcribed_processed_pseudogene']
df_biotype = df_biotype[df_biotype['biotype'].isin(accepted)]

# From df_parent, keeping only pseudogenes that have been processed
df_parent = df_parent[df_parent['ID'].isin(df_biotype['transcript_id'].to_list())]

# Couting pseudogenes per gene
pseudo_count = df_parent.groupby('Parent gene').count()['ID'].sort_values()
pseudo_count.index = pseudo_count.index.str.split('.').str[0]
pseudo_dict = pseudo_count.to_dict()

# List of genes present in dataset and not already pseudogenes
genes = df_gene_biotype[
    ~df_gene_biotype['gene_biotype'].str.contains('pseudo')
]['gene_id'].tolist()

# Keeping only genes present in the dataset
ens_genes = [
    gene for gene in data.index.get_level_values('ensembl_id').to_list()
    if gene in genes
]
pseudo_count = pd.DataFrame(index=ens_genes)
pseudo_count['pseudogene'] = pseudo_count.index.map(pseudo_dict)
pseudo_count.fillna(0, inplace=True)

# Plotting
plt.figure(figsize=mm2inch((40,50)))
plt.boxplot(
    [
        pseudo_count[pseudo_count.index.isin(lgenes)]['pseudogene'],
        pseudo_count[pseudo_count.index.isin(lgenes2)]['pseudogene']
    ],
    widths=.8
)

plt.ylim([-0.1, 160])
plt.ylabel('Number of processed pseudogene')
plt.yscale('function', functions=(np.arcsinh, np.sinh))
plt.yticks([0, 1, 3, 7, 15, 30, 75, 150])
plt.tight_layout()
plt.savefig(snakemake.output.plot)
