import numpy as np
import pandas as pd
import snakemake.shell as shell

# Loading parameters
columns = [
    'seq', 'source', 'feature', 'start', 'end',
    'score', 'strand', 'frame', 'attributes'
]

def get_gene(gtf_file, gene):
    # Get GTF information for gene
    gtf_gene = shell(
        "cat {gtf_file} | grep '\"{gene}\"' ",
        iterable=True
    )

    # Format GTF file
    gtf_lines = [line.split('\t') for line in gtf_gene]
    df_gene = pd.DataFrame(gtf_lines, columns=columns)

    # Adding columns
    df_gene['transcript_name'] = df_gene['attributes'].str.extract("transcript_name \"(.*?)\"")

    # Extracting features for the gene
    gene_indx = df_gene['feature'] == 'gene'
    if np.sum(gene_indx) == 1:
        gene_start = int(df_gene.loc[gene_indx, 'start'][0])
        gene_end = int(df_gene.loc[gene_indx, 'end'][0])
        gene_strand = df_gene.loc[gene_indx, 'strand'][0]
        gene_seq = df_gene.loc[gene_indx, 'seq'][0]
    else:
        print('More than one gene in gtf using symbol ' + snakemake.wildcards.gene)
        raise AssertionError

    return df_gene, gene_start, gene_end, gene_strand, gene_seq
