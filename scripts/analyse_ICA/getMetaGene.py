import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from get_gene_gtf import get_gene

# Params
min_length = 160
box_length = 80
extending = 10

box_length = 120
min_length = 2*box_length
extending = 15

primes = [3, 5]
dist = dict()
for prime in primes:
    dist[prime] = {
        "HISAT2": [],
        "STAR": [],
        "tophat2": []
    }

for gene_file in snakemake.input.genes:
    gene_dist = dict()
    for prime in primes:
        gene_dist[prime] = {
            "HISAT2": [],
            "STAR": [],
            "tophat2": []
        }
    gene = gene_file.split('/')[-1].split('.')[0]
    df_gene, gene_start, gene_end, gene_strand, gene_seq = get_gene(
        snakemake.input.gtf, gene
    )

    # Loading data
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
        dist_3 = None

        # Exon parameters
        start = int(exon.start) - gene_start
        end = int(exon.end) - gene_start
        exon_len = end-start
        if exon_len > min_length:

            # For positive strand
            if gene_strand == '+':
                # Not considering the 5' end of the first exon
                if int(exon.exon_number) != 1:
                    dist_5 = data.iloc[start-extending-1:start+box_length]
                # Not considering the 3' end of the last exon
                if int(exon.exon_number) != transc_count[exon.transcript_name]:
                    dist_3 = data.iloc[end-box_length:end+extending+1]

            # For negative strand
            elif gene_strand == '-':
                if int(exon.exon_number) != transc_count[exon.transcript_name]:
                    dist_5 = data.iloc[end-box_length-1:end+extending][::-1]
                if int(exon.exon_number) != 1:
                    dist_3 = data.iloc[start-extending:start+box_length+1][::-1]

            for _dist, prime in [(dist_5, 5), (dist_3, 3)]:
                if _dist is not None:
                    max = np.max([np.max(_dist[tool]) for tool in _dist.keys()])
                    for tool in gene_dist[prime].keys():
                        _dist[tool] /= max
                        gene_dist[prime][tool].append(_dist[tool].tolist())

    # Getting the mean profile for each gene
    for prime in primes:
        for tool in gene_dist[prime].keys():
            data = pd.DataFrame(gene_dist[prime][tool])
            dist[prime][tool].append(data.mean(axis=0).tolist())

data = list()
cols = list()
for prime in dist.keys():
    for tool in dist[prime].keys():
        for indx, genes in enumerate(dist[prime][tool]):
            cols.append((prime, tool, indx))
            data.append(genes)

data = pd.DataFrame(data, index=pd.MultiIndex.from_tuples(cols, names=['prime', 'tool', 'gene'])).transpose()
mean = data.groupby(level=['prime', 'tool'], axis=1).mean()
std = data.groupby(level=['prime', 'tool'], axis=1).std()

# To file
mean.to_csv(snakemake.output.mean, sep='\t')
std.to_csv(snakemake.output.std, sep='\t')
