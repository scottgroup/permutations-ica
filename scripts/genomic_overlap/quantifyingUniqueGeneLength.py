import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
import seaborn as sns

from matplotlib.patches import Rectangle

from plotting_utils import mm2inch, rcParams, reading_file, colors
for k, v in rcParams.items():
    plt.rcParams[k] = v


# Loading HGNC data
hgnc_df = pd.read_csv(snakemake.input.HGNC, sep='\t')
hgnc_df['ncbi_id'] = "GeneID:" + hgnc_df['ncbi_id'].astype('str').str.split('.').str[0]

# Creating translation dictionary
ref2hgnc = hgnc_df.set_index('ncbi_id')['HGNC_ID'].to_dict()
ens2hgnc = hgnc_df.set_index('ensembl_id')['HGNC_ID'].to_dict()

# Loading data
annotations = ['ensembl92', 'ensembl98', 'refseq']

# Loading gene list
gene_list_fname = snakemake.params.gene_list_dir + 'comp_{comp}.txt'.format(
    comp=snakemake.wildcards.comp
)
up, down = reading_file(gene_list_fname)
up = [ens2hgnc[gene] for gene in up]
down = [ens2hgnc[gene] for gene in down]

# Checking out annotations
fig, axes = plt.subplots(nrows=1, ncols=2, figsize=mm2inch(60, 50))
j = 0
for ori, ori_str in zip([down, up], ['down', 'up']):
    len_annot = {annot:0 for annot in annotations}
    res = list()
    errors = list()
    for i, annotation in enumerate(snakemake.input.annotations):
        df = pd.read_csv(annotation, sep='\t', names=['gene1', 'gene2', 'score1', 'score2'])
        annot1, annot2 = annotation.split('_')[-3], annotation.split('_')[-2]

        # Translating everything to HGNC_ID
        for gene in ['gene1', 'gene2']:
            if 'ENSG' in df[gene].tolist()[0]:
                df[gene] = df[gene].map(ens2hgnc)
            elif 'GeneID' in df[gene].tolist()[0]:
                df[gene] = df[gene].map(ref2hgnc)
        df.dropna(axis=0, inplace=True)

        # Keeping only scores for the same gene
        df = df[df['gene1'] == df['gene2']]

        # Keeping only genes from orientation
        scoring_df = df[df['gene1'].isin(ori)].set_index('gene1')
        scoring_df = pd.DataFrame(scoring_df, index=ori)[['score1', 'score2']]
        scoring_df.fillna(value=0, inplace=True)

        score1, score2 = scoring_df.mean(axis=0).tolist()
        std1, std2 = scoring_df.std(axis=0).tolist()

        # Vectors are
        # [refseq, Rvs98, 98, 98vs92, Rv92, 92]
        data = [0, 0, 0, 0, 0, 0]
        error = [0, 0, 0, 0, 0, 0]
        annot_indx = {"refseq": 0, "ensembl98": 2, "ensembl92": 5}

        _len = 1
        lens = [_len/score1 - _len, _len, _len/score2 - _len]
        lens = [_len/sum(lens) for _len in lens]

        for i, annot in enumerate(annotations):
            if annot1 == annot:
                data[annot_indx[annot]] = lens[0]
                error[annot_indx[annot]-1] = std1/score1 * lens[0]
            elif annot2 == annot:
                data[annot_indx[annot]] = lens[2]
                error[annot_indx[annot]] = std2/score2 * lens[2]

            if "refseq" in [annot1, annot2]:
                if "ensembl98" in [annot1, annot2]:
                    data[1] = lens[1]
                elif "ensembl92" in [annot1, annot2]:
                    data[4] = lens[1]
            else:
                data[3] = lens[1]

        res.append(data)
        errors.append(error)

    # Plotting!
    bottom = [0, 0, 0]
    for i in range(len(res[0])):
        data = [datum[i] for datum in res]
        err = [_err[i] for _err in errors]
        color='k'
        if i == 0:
            color = colors['refseq']
        elif i == 2:
            color = colors['ensembl98']
        elif i == 5:
            color = colors['ensembl92']
        else:
            color = 'gray'

        axes[j].bar(
            [0, 1, 2], data, color=color, bottom=bottom, yerr=err
        )
        bottom = list(np.array(bottom) + np.array(data))

    j += 1
plt.tight_layout()
plt.savefig(snakemake.output.plot)
