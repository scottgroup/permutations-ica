import os
import pandas as pd
from snakemake import shell

# Importing parameters
annotation = snakemake.wildcards.annotation
gene = snakemake.wildcards.gene
gtf = snakemake.input.gtf
exon_bed = snakemake.input.exon_bed

# Creating temp bed file
gene_bed = snakemake.output.bed
gene_bed_merged = snakemake.output.merged_bed

# Transforming gene for annotation specific format
if 'ENSG' in gene:
    gene_id = 'gene_id "{gene}"'.format(gene=gene)
else:
    gene_id = 'db_xref "{gene}"'.format(gene=gene)

# Getting new gene_id
gene_id = shell(
    "zcat {gtf} | grep '{gene_id}' ",
    iterable=True
)
gene_id = list(gene_id)[0].split('\t')[-1].split(';')[0]

# Verify if Refseq is not doing its thing where a gene does not have exon
check = shell("zcat {gtf} | grep '{gene_id}'", iterable=True)
if len(list(check)) == 1:
    zcat_command =  "zcat {gtf} | grep '{gene_id}' | "
else:
    zcat_command = "zcat {gtf} | grep '{gene_id}' | grep '\texon\t' | "

# Reading through the gtf file and write merged gene bed
shell(
    zcat_command +
    "awk '{{print $1,$4-1,$5,\"{gene}\",1,$7}}' OFS=\"\t\" FS=\"[\t;]\" | "
    "sort -k1,1 -k2,2n > {gene_bed}"
)
shell(
    "bedtools merge -i {gene_bed} -c 4,5,6 -o distinct,distinct,distinct "
    "> {gene_bed_merged}"
)

# Get the overlap between the merged gene bed and the exons bed
other_genes = shell(
    "cat {exon_bed} | grep -v {gene} | "
    "bedtools intersect -a {gene_bed_merged} -b stdin | sort -k1,1 -k2,2n | bedtools merge",
    iterable=True
)

columns = ['chr', 'start', 'end', 'id', 'score', 'strand']
df = pd.read_csv(gene_bed_merged, sep='\t', names=columns)
o_df = pd.DataFrame(
    [other_gene.split('\t') for other_gene in other_genes],
    columns=columns[:3]
)

r_gene = abs(sum(df['end'].astype('int') - df['start'].astype('int')))
r_other = abs(sum(o_df['end'].astype('int') - o_df['start'].astype('int')))
res = r_other/r_gene

with open(snakemake.output.overlap, 'w') as f:
    f.write('\t'.join([annotation, gene, str(res)]) + '\n')
