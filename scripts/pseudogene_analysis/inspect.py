import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats

# Using data from http://pseudogene.org/psicube/index.html

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
    snakemake.input.data,
    sep='\t', usecols=[0,1,2,3], skiprows=6
)

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

# Keeping only genes present in the dataset
ens_genes = data['ensembl_id'].to_list()
pseudo_count = pseudo_count[pseudo_count.index.isin(ens_genes)]

print(pseudo_count)


# Loading other
pseudo_blat = pd.read_csv(snakemake.input.stuff, sep='\t')
blat_count = pseudo_blat.groupby('subjectGene').count()['query'].sort_values()
blat_count = blat_count[blat_count.index.isin(ens_genes)]

print(blat_count)
print(pseudo_count[pseudo_count.index.isin(up)])
print(blat_count[blat_count.index.isin(up)])

print('There is ', len(up), ' up genes')

# stuff = (blat_count - pseudo_count)
# print(stuff.sort_values())

# Plotting?
plt.boxplot(
    [
        pseudo_count[~pseudo_count.index.isin(up)],
        pseudo_count[pseudo_count.index.isin(up)],
        blat_count[~blat_count.index.isin(up)],
        blat_count[blat_count.index.isin(up)],    ],
    widths=.8
)
plt.ylabel('Number of processed pseudogene')
plt.show()

qwe = stats.ttest_ind(
    pseudo_count[~pseudo_count.index.isin(up)],
    pseudo_count[pseudo_count.index.isin(up)],
    equal_var=False
)
print(qwe)

qwe = stats.ttest_ind(
    blat_count[~blat_count.index.isin(up)],
    blat_count[blat_count.index.isin(up)],
    equal_var=False
)
print(qwe)
