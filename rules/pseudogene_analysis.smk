
def get_counts(wildcards):
    """ """
    fname = config['ICA_datasets'][wildcards.dataset]['params']['counts']
    return rna_seq_cartesian_product("results/cartesian_product/{dataset}.tsv".format(dataset=fname))


rule download_pseudogene_parents:
    output:
        file = "data/references/pseudogene_parents.txt"
    params:
        link = lambda w: "{pseudogenes_parents}".format(**config['url'])
    shell:
        "wget --quiet -O {output.file} {params.link}"


rule download_pseudogene_biotype:
    output:
        file = "data/references/pseudogene_biotype.txt"
    params:
        link = lambda w: "{pseudogenes_biotype}".format(**config['url'])
    shell:
        "wget --quiet -O {output.file} {params.link}"


rule get_gene_biotype:
    input:
        gtf = rna_seq_cartesian_product("data/references/ensembl98.gtf")
    output:
        biotype = "data/references/gene_biotypes.tsv"
    conda:
        "../envs/ICA_python.yaml"
    script:
        "../scripts/pseudogene_analysis/get_gene_biotype.py"


rule inspect_pseudogenes_parents:
    input:
        parent = rules.download_pseudogene_parents.output.file,
        biotype = rules.download_pseudogene_biotype.output.file,
        gene_list_tkn = "results/ICA/{ICAmethod}/{dataset}/{ICA_run}/gene_list/comp_sigma{sigma}/.tkn",
        data = "subworkflows/rna_seq_cartesian_product/results/cartesian_product/tissues_NaN.tsv", # get_counts,
        gene_biotype = rules.get_gene_biotype.output.biotype
    output:
        plot = "results/ICA/{ICAmethod}/{dataset}/{ICA_run}/sigma_{sigma}/pseudogene_comp{comp}.svg"
    params:
        gene_list = "results/ICA/{ICAmethod}/{dataset}/{ICA_run}/gene_list/comp_sigma{sigma}/comp_{comp}.txt"
    conda:
        "../envs/ICA_python.yaml"
    script:
        "../scripts/pseudogene_analysis/inspect.py"
