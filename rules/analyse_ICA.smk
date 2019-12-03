
def getGeneBEDs(wildcards):

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

    BEDs = list()
    path = "results/rnaseq/geneCoverage/{gene}.tsv"

    for gene in genes:
        BEDs.append(rna_seq_cartesian_product(path.format(gene=gene)))

    return BEDs


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


rule plotMetaGene:
    input:
        genes = getGeneBEDs
    output:
        plot = "results/{ICA_path}/comp_{comp}_{side}_sigma{sigma}.svg"
    conda:
        "../envs/ICA_python.yaml"
    script:
        "../scripts/analyse_ICA/plotMetaGene.py"
