import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats

from plotting_utils import mm2inch, rcParams
for k, v in rcParams.items():
    plt.rcParams[k] = v

# Loading HGNC
hgnc = pd.read_csv(snakemake.input.HGNC, sep='\t')
for c in hgnc.columns:
    hgnc[c] = hgnc[c].astype('str')
    hgnc[c].apply(str)
hgnc['ncbi_id'] = 'GeneID:' + hgnc['ncbi_id'].str.split('.').str[0]

for overlap_path in snakemake.input.overlaps:
    annotation = overlap_path.split('/')[-1].split('.')[0]
    if 'ensembl' in annotation:
        col = 'ensembl_id'
    else:
        col = 'ncbi_id'

    overlaps = pd.read_csv(overlap_path, sep='\t', index_col=0)
    over_dict = overlaps['score'].to_dict()
    hgnc[annotation] = hgnc[col].map(over_dict)

hgnc = hgnc.dropna()
hgnc.drop(labels=['symbol', 'ncbi_id', 'ensembl_id'], axis=1, inplace=True)
hgnc.set_index('HGNC_ID', inplace=True)

# Removing 0
data = hgnc.melt()
data = data[data['value'] > 0]

# Plotting violin
fig = plt.figure(figsize=mm2inch((80, 45)))
sns.violinplot(y='variable', x='value', data=data, cut=0, bw=0.25, scale='count')
plt.xlabel('% of exon overlapped')
plt.ylabel('')
plt.xlim([0,1])
plt.tight_layout()
plt.savefig(snakemake.output.plot)
