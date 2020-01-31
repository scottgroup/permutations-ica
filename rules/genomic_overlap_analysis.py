import pandas as pd
import os

annotations = ["ensembl92", "ensembl98", "refseq"]
count_file = "subworkflows/rna_seq_cartesian_product/results/cartesian_product/tissues_NaN.tsv"
HGNC_path = "subworkflows/rna_seq_cartesian_product/data/references/hgnc.txt"
path_adj = "results/genomic_overlap/adjacent_genes/{gene}_{annotation}.tsv"
bed_path = "results/genomic_overlap/gene_bed_merged/{gene}_{annotation}.bed"
rt_overlap_path = "results/genomic_overlap/readthrough/{HGNCgene}.tsv"
overlap_path = "results/genomic_overlap/gene_overlap/{gene}_{annotation}.tsv"


def getAdjacentGenes(path, wildcards):
    """
        For a specific gene, specified with an HGNC ID :

        Ask for path_adj of the gene, a file which contains all genes sharing a
        genomic coordinate with the gene of interest.

        If path_adj exist, require the (path=path_adj or bed_path) of the genes
        inside that path_adj of the HGNC gene.

        If path == 'adj' and adjacent gene of the gene of interest have been
        found, find their adjacents genes. The objectif is to detect overlapping
        read-throughts by looking at how many other genes they overlap.
    """
    if path == 'adj':
        fpath = path_adj
    elif path == 'bed':
        fpath = bed_path

    # Loading HGNC and relevant dictionaries
    HGNC = wildcards['HGNCgene']
    _, HGNC2ens, HGNC2ref = getHGNCdict(HGNC_path)
    ens_gene = HGNC2ens[HGNC]
    ref_gene = HGNC2ref[HGNC]
    gene2HGNC = {
        **{val:key for key, val in HGNC2ens.items()},
        **{val:key for key, val in HGNC2ref.items()}
    }

    # Loading list of adjacent genes
    files = list()
    for annotation in annotations:
        if 'ensembl' in annotation:
            gene = ens_gene
        else:
            gene = ref_gene

        file = path_adj.format(gene=gene,annotation=annotation)
        if os.path.exists(file):
            with open(file, 'r') as f:
                for line in f.readlines():
                    _gene, _annot = line.strip().split('\t')
                    _file = fpath.format(gene=_gene,annotation=_annot)
                    files.append(_file)

                    if os.path.exists(_file) and path == 'adj':
                        with open(_file, 'r') as f:
                            for line in f.readlines():
                                _gene, _annot = line.strip().split('\t')
                                _file = bed_path.format(gene=_gene,annotation=_annot)
                                files.append(_file)
        else:
            files.append(file)

    return list(set(files))


def getGeneOverlaps(wildcards):
    """  From a gene list query all gene overlap files for each annotation """
    path = overlap_path
    comp_path = "results/{ICA_path}/sigma_{sigma}/gene_list/comp_{comp}.txt"
    genes = list()
    with open(comp_path.format(**wildcards), 'r') as f:
        for line in f.readlines():
            line = line.strip()
            if line[0] != '>':
                genes.append(line)

    ens2ref, _, _ = getHGNCdict(HGNC_path)

    files = list()
    for gene in genes:
        for annotation in annotations:
            if annotation is "refseq":
                gene = ens2ref[gene]
            files.append(path.format(gene=gene, annotation=annotation))
    return files


def getHGNCdict(HGNC_path):
    """ Create dictionaries linking between gene IDs """
    # Loading HGNC data
    hgnc_df = pd.read_csv(HGNC_path, sep='\t')
    hgnc_df['ncbi_id'] = "GeneID:" + hgnc_df['ncbi_id'].astype('str').str.split('.').str[0]

    # Creating translation dictionary
    ens2ref = hgnc_df.set_index('ensembl_id')['ncbi_id'].to_dict()
    HGNC2ref = hgnc_df.set_index('HGNC_ID')['ncbi_id'].to_dict()
    HGNC2ens = hgnc_df.set_index('HGNC_ID')['ensembl_id'].to_dict()

    return ens2ref, HGNC2ens, HGNC2ref
