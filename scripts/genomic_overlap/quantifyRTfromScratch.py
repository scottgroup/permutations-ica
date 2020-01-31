import pandas as pd
import itertools

column_types = {
    'seq': 'object',
    'feature': 'object',
    'start': 'int64',
    'end': 'int64',
    'strand': 'object',
    'attributes': 'object',
}

def process_genes(genes, hit_min=1):
    results = list()
    for gene in genes:
        if str(gene)[0] not in ['E', 'G']:
            gene = 'GeneID:' + str(gene)
        gene_gtf = gtf[gtf['gene_id'] == gene]
        hits = dict()

        chrs = gene_gtf['seq'].unique()
        strands = gene_gtf['strand'].unique()
        sub_gtf = gtf[gtf['seq'].isin(chrs) & gtf['strand'].isin(strands)]

        for idx, row in gene_gtf.iterrows():
            seq, start, end, strand, gene_id = row
            _h = sub_gtf[(sub_gtf['start'] == start) &(sub_gtf['end'] == end)]

            if len(_h) > 1:
                _genes_hit = [g for g in _h['gene_id'].to_list() if g != gene]
                for _gene_hit in _genes_hit:
                    if _gene_hit not in hits:
                        hits[_gene_hit] = 0
                    hits[_gene_hit] += 1

        for k, v in hits.items():
            if v >= hit_min:
                results.append(
                    [gene, k, v]
                )

    return pd.DataFrame(results, columns=['gene', 'overlap', 'N_exon'])


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
gtf = pd.read_csv(
    snakemake.input.gtf, sep='\t', skiprows=[0,1,2,3,4], usecols=[0,2,3,4,6,8],
    names=[
        'seq', 'source', 'feature', 'start', 'end',
        'score', 'strand', 'frame', 'attributes'
    ]
)
gtf.dropna(axis=0, inplace=True)
gtf = gtf.astype(column_types)

gtf['gene_id'] = gtf['attributes'].str.extract(extract)
gtf['gene_id'] = gtf['gene_id'].fillna(method='ffill')
if snakemake.wildcards.annotation == "refseq":
    gtf['gene_id'] = 'GeneID:' + gtf['gene_id']
gtf.drop(columns='attributes', inplace=True)

# Keeping only exon
gtf = gtf[gtf['feature'] == 'exon']
gtf.drop(columns='feature', inplace=True)

# Removing duplicates
gtf.drop_duplicates(inplace=True)

# For all genes
overlaps = process_genes(genes)

# Finding genes that need completion
to_quant = [g for g in overlaps['overlap'] if g not in overlaps['gene']]
_overlaps = process_genes(to_quant)
overlaps = pd.concat([overlaps, _overlaps])

# Identifying potential read-throughs
n_overlap = overlaps.groupby('gene').count().reset_index()
read_through = n_overlap[n_overlap['overlap'] > 1]

# Need to confirm that they use different exons
possible_RT = read_through['gene'].unique()
RTs = process_genes(possible_RT)
real_RTs = list()
for gene in RTs['gene'].unique():
    hit = False
    testOverlap = RTs[RTs['gene'] == gene]['overlap'].tolist()
    for _is_overlap in testOverlap:
        resp = process_genes([_is_overlap])['overlap'].tolist()
        for to_overlap in [g for g in testOverlap if g != _is_overlap]:
            if to_overlap not in resp:
                hit = True
    if hit:
        real_RTs.append(gene)

RTs = RTs[RTs['gene'].isin(real_RTs)]
RTs.to_csv(snakemake.output.RT_list, sep='\t', index=False)
