import pandas as pd
import numpy as np

column_types = {
    'seq': 'object',
    'feature': 'object',
    'start': 'int64',
    'end': 'int64',
    'strand': 'object',
    'attributes': 'object',
}

# Checking wilcards
if snakemake.wildcards.annotation == "refseq":
    gene_id = "ncbi_id"
    extract = "db_xref \"GeneID:(.*?)\""
else:
    gene_id = "ensembl_id"
    extract = "gene_id \"(.*?)\""

# Loading genes
genes = pd.read_csv(
    snakemake.input.count_file, sep='\t', usecols=[0,2,3], skiprows=[0,1,2,3,4,5]
)
genes = genes[gene_id].to_list()

# Loading annotation
gtf = pd.read_csv(
    snakemake.input.gtf, sep='\t', skiprows=[0,1,2,3,4], usecols=[0,2,3,4,6,8],
    names=[
        'seq', 'source', 'feature', 'start', 'end',
        'score', 'strand', 'frame', 'attributes'
    ]
)
gtf.dropna(axis=0, inplace=True)
gtf = gtf.astype(column_types)

# Correction for open/closed end
gtf['end'] = gtf['end'] + 1

gtf['gene_id'] = gtf['attributes'].str.extract(extract)
gtf['gene_id'] = gtf['gene_id'].fillna(method='ffill')
if snakemake.wildcards.annotation == "refseq":
    gtf['gene_id'] = 'GeneID:' + gtf['gene_id']
gtf.drop(columns='attributes', inplace=True)

# Keeping only exon
gtf_gene = gtf[gtf['feature'] == 'gene']
gtf = gtf[gtf['feature'] == 'exon']
gtf.drop(columns='feature', inplace=True)

# Removing duplicates
gtf.drop_duplicates(inplace=True)

# Creating a dict for each seq
chr_dict = dict()
for chr in gtf['seq'].unique():
    chr_dict[chr] = gtf[gtf['seq'] == chr]

# Crawling the genes
overlaps = list()
for gene in genes:
    if snakemake.wildcards.annotation == "refseq":
        gene = "GeneID:" + str(gene)
    _gene = gtf_gene[gtf_gene['gene_id'] == gene]

    if len(_gene) >= 1:
        seq, feature, start, end, strand, gene_id  = _gene.values[0]

        # Creating gene vector
        gene_loc = np.zeros(end-start)
        other_loc = np.zeros(end-start)

        # Extracting exon_gtf from the gene
        exon_gtf = chr_dict[seq][chr_dict[seq]['gene_id'] == gene]
        for idx, exon in exon_gtf.iterrows():
            gene_loc[exon['start']-start:exon['end']-start] = 1

        # Extracting exon_gtf for the other genes
        other_gtf = chr_dict[seq][
            (chr_dict[seq]['start'] < end) &
            (chr_dict[seq]['end'] > start) &
            (chr_dict[seq]['gene_id'] != gene)
        ]
        other_gtf.loc[other_gtf['end'] > end, 'end'] = end
        other_gtf.loc[other_gtf['start'] < start, 'start'] = start

        for idx, exon in other_gtf.iterrows():
            other_loc[exon['start']-start:exon['end']-start] = 1

        intersection = gene_loc * other_loc
        score = np.sum(intersection)/np.sum(gene_loc)
        overlaps.append((gene, score))

df = pd.DataFrame(overlaps, columns=['gene', 'score'])
df.to_csv(snakemake.output.overlaps, sep='\t', index=False)
