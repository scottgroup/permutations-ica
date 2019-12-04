from functools import partial


def getGeneFiles(path, wildcards):
    """ """
    gene_list = "results/{ICA_path}/gene_list/comp_sigma{sigma}/comp_{comp}.txt"
    genes = list()
    with open(gene_list.format(**wildcards), 'r') as f:
        for line in f.readlines():
            if line.strip() == ">Positive genes":
                up = True
            elif line.strip() == ">Negative genes":
                up = False
            elif wildcards.side == 'up' and up:
                genes.append(line.strip())
            elif wildcards.side == 'down' and not up:
                genes.append(line.strip())

    files = list()
    for gene in genes:
        files.append(rna_seq_cartesian_product(path.format(gene=gene)))
    return files


rule GO_analysis:
    input:
        gene_list = "results/{ICA_path}/gene_list/comp_sigma{sigma}/.tkn",
        components_mean = "results/{ICA_path}/components_mean.tsv"
    output:
        plot = "results/{ICA_path}/gene_list/comp_sigma{sigma}/GO_comp{comp}_{ori}.svg",
    params:
        path = "results/{ICA_path}/gene_list/comp_sigma{sigma}/comp_{comp}.txt"
    conda:
        "../envs/gprofiler.yaml"
    script:
        "../scripts/analyse_ICA/GO_analysis.py"


rule getMetaGene:
    input:
        genes = partial(getGeneFiles, "results/rnaseq/geneCoverage/{gene}.tsv"),
        gtf = rna_seq_cartesian_product("data/references/ensembl98.gtf")
    output:
        mean = "results/{ICA_path}/metaGene/comp_{comp}_{side}_sigma{sigma}_mean.tsv",
        std = "results/{ICA_path}/metaGene/comp_{comp}_{side}_sigma{sigma}_std.tsv"
    conda:
        "../envs/ICA_python.yaml"
    script:
        "../scripts/analyse_ICA/getMetaGene.py"


rule plotMetaGene:
    input:
        mean = rules.getMetaGene.output.mean,
        std = rules.getMetaGene.output.std
    output:
        plot = "results/{ICA_path}/metaGene/comp_{comp}_{side}_sigma{sigma}.svg"
    conda:
        "../envs/ICA_python.yaml"
    script:
        "../scripts/analyse_ICA/plotMetaGene.py"


rule getOtherChr:
    input:
        gene_BAM_tkn = partial(getGeneFiles, "results/geneCoverage/{gene}.tkn"),
    output:
        file = "results/{ICA_path}/metaGene/comp_{comp}_{side}_sigma{sigma}_chr"
