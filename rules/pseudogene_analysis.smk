
def get_counts(wildcards):
    """ Returns path of the count file """
    fname = config['ICA_models'][wildcards.ICAmodel]['params']['counts']
    return rna_seq_cartesian_product("results/cartesian_product/{dataset}.tsv".format(dataset=fname))


rule download_pseudogene_parents:
    """ Download pseudogene data from the PsiCube project """
    output:
        file = "data/references/pseudogene_parents.txt"
    params:
        link = lambda w: "{pseudogenes_parents}".format(**config['url'])
    shell:
        "wget --quiet -O {output.file} {params.link}"


rule download_pseudogene_biotype:
    """ Download pseudogene parents biotype data from the PsiCube project """
    output:
        file = "data/references/pseudogene_biotype.txt"
    params:
        link = lambda w: "{pseudogenes_biotype}".format(**config['url'])
    shell:
        "wget --quiet -O {output.file} {params.link}"


rule get_gene_biotype:
    """ Extracts biotype from Ensembl98 genes """
    input:
        gtf = rna_seq_cartesian_product("data/references/ensembl98.gtf")
    output:
        biotype = "data/references/gene_biotypes.tsv"
    conda:
        "../envs/ICA_python.yaml"
    script:
        "../scripts/pseudogene_analysis/get_gene_biotype.py"


rule compare_pseudogenes_parents:
    """
        Plots a distribution of the number of pseudogenes for each gene in two
        compared components.
    """
    input:
        parent = rules.download_pseudogene_parents.output.file,
        biotype = rules.download_pseudogene_biotype.output.file,
        gene_list_tkn = "results/ICA/{ICAmethod}/{ICAmodel}/{ICA_run}/sigma_{sigma}/gene_list/.tkn",
        data = get_counts,
        gene_biotype = rules.get_gene_biotype.output.biotype
    output:
        plot = "results/ICA/{ICAmethod}/{ICAmodel}/{ICA_run}/sigma_{sigma}/pseudogene_comp{comp}_vs_comp{comp2}.svg"
    params:
        gene_list = "results/ICA/{ICAmethod}/{ICAmodel}/{ICA_run}/sigma_{sigma}/gene_list/comp_{comp}.txt",
        gene_list2 = "results/ICA/{ICAmethod}/{ICAmodel}/{ICA_run}/sigma_{sigma}/gene_list/comp_{comp2}.txt"
    conda:
        "../envs/ICA_python.yaml"
    script:
        "../scripts/pseudogene_analysis/compare_pseudogenes_parents.py"
