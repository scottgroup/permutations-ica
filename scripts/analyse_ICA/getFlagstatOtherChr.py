import pandas as pd
import snakemake.shell as shell

def get_bam_variables(fname, config):
    """ """
    bam_var = dict()

    for tool, vars in config['tools'].items():
        for var in vars:
            if var in fname.split('/'):
                bam_var[tool] = var

    for tissue, datasets in config['datasets'].items():
        if tissue in fname:
            bam_var['tissue'] = tissue
            for dataset in datasets:
                if dataset in fname:
                    bam_var['dataset'] = dataset

    return bam_var

# Get all genes
genes = [
    gene.split('/')[-1].split('.')[0]
    for gene in snakemake.input.gene_BAM_tkn
]

# Get all folders
folders = snakemake.input.folders

header = list()
genes_quant = list()

for folder in folders:
    gene_quant = list()
    header.append(get_bam_variables(folder, snakemake.config))
    for gene in genes:
        flagstats = list(
            shell("samtools flagstat {folder}/{gene}.bam", iterable=True)
        )
        paired = int(flagstats[5].split()[0])
        other_chr = int(flagstats[11].split()[0])
        if paired == 0:
            perc = 0
        else:
            perc = other_chr/paired*100
        gene_quant.append(perc)

    genes_quant.append(gene_quant)

df_header = pd.DataFrame(header)
df = pd.DataFrame(genes_quant, index=pd.MultiIndex.from_frame(df_header), columns=genes).T
df = df.groupby(level=['aligner', 'dataset', 'tissue', 'trimmer'], axis=1).mean()
df.to_csv(snakemake.output.file, sep='\t')
