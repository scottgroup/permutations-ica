import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats

plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams['font.size'] = 6
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial']


def mm2inch(*tupl):
    """
    https://stackoverflow.com/questions/14708695/specify-figure-size-in-centimeter-in-matplotlib
    """
    inch = 25.4
    if isinstance(tupl[0], tuple):
        return tuple(i/inch for i in tupl[0])
    else:
        return tuple(i/inch for i in tupl)


def reading_file(path):
    up, down = list(), list()
    with open(path, 'r') as f:
        for line in f.readlines():
            line = line.strip()
            if line == '>Positive genes':
                reading_up = True
            elif line == '>Negative genes':
                reading_up = False
            else:
                if reading_up:
                    up.append(line)
                else:
                    down.append(line)
    return up, down

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
up, down = reading_file(snakemake.params.gene_list)

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
plt.figure(figsize=mm2inch((40,65)))
plt.boxplot(
    [
        pseudo_count[pseudo_count.index.isin(up)]['pseudogene'],
        pseudo_count[~pseudo_count.index.isin(up)]['pseudogene']
    ],
    widths=.8
)
plt.ylim([0,159])
plt.ylabel('Number of processed pseudogene')
plt.tight_layout()
plt.savefig(snakemake.output.plot)

stat = stats.ttest_ind(
    pseudo_count[pseudo_count.index.isin(up)],
    pseudo_count[~pseudo_count.index.isin(up)],
    equal_var=False
)
print(stat)
