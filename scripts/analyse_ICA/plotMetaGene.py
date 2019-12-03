import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from get_gene_gtf import get_gene


# Params
min_length = 160
box_length = 80


dist = {
    "HISAT2": [],
    "STAR": [],
    "tophat2": []
}

gene_dist = {
    "HISAT2": [],
    "STAR": [],
    "tophat2": []
}

for gene_file in snakemake.input.genes:
    gene = gene_file.split('/')[-1].split('.')[0]
    print(gene)

    df_gene, gene_start, gene_end, gene_strand, gene_seq = get_gene(
        snakemake.input.gtf, gene
    )
    print(gene_strand)

    data = pd.read_csv(gene_file, sep='\t', index_col=0, header=[0, 1, 2, 3, 4])
    data = data.groupby(level=['aligner', 'dataset', 'tissue', 'trimmer'], axis=1).mean()
    data = data.groupby(level='aligner', axis=1).mean()

    # For all exon
    df_exon = df_gene[df_gene['feature'] == 'exon']
    transcripts = df_exon['transcript_name'].unique().tolist()
    transc_count = {
        transcript: len(df_exon[df_exon['transcript_name'] == transcript])
        for transcript in transcripts
    }

    for indx, exon in df_exon.iterrows():
        # Reseting
        dist_5 = None
        # print(exon.start, exon.end)
        start = int(exon.start) - gene_start
        end = int(exon.end) - gene_start

        exon_len = end-start
        if exon_len > min_length:

            # For positive strand
            if gene_strand == '+':
                # Not considering the 5' end of the first exon
                if int(exon.exon_number) != 1:
                    dist_5 = data.iloc[start-5:start+box_length]

            # For negative strand
            elif gene_strand == '-':
                pass

            if dist_5 is not None:
                max = np.max([np.max(dist_5[tool]) for tool in dist_5.keys()])
                for tool in gene_dist.keys():
                    dist_5[tool] /= max
                    gene_dist[tool].append(dist_5[tool].tolist())

    # Getting the mean profile for each gene
    for tool in gene_dist.keys():
        data = pd.DataFrame(gene_dist[tool])
        dist[tool].append(data.mean(axis=0).tolist())


for tool in dist.keys():
    mean = np.array(pd.DataFrame(dist[tool]).mean(axis=0).tolist())
    std = np.array(pd.DataFrame(dist[tool]).std(axis=0).tolist())

    x = np.arange(len(mean))

    plt.plot(x, mean, linewidth=2)
    plt.fill_between(
    x = x,
    y1 = mean + std,
    y2 = mean - std,
    alpha=0.3
)

plt.savefig(snakemake.output.plot)
