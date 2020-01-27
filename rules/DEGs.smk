from DEGs import get_DEG_results


def get_counts(wildcards):
    """ """
    fname = config['ICA_models'][wildcards.ICAmodel]['params']['counts']
    return rna_seq_cartesian_product("results/cartesian_product/{dataset}.tsv".format(dataset=fname))


rule prepare_data_for_DESeq2:
    """ Format the counts and samples matrices for DESeq2 """
    input:
        dataset = get_counts,
    output:
        counts = "results/DESeq2/{ICAmodel}/counts.tsv",
        samples = "results/DESeq2/{ICAmodel}/samples.tsv"
    conda:
        "../envs/ICA_python.yaml"
    script:
        "../scripts/DEGs/1_prepare_data_for_DESeq2.py"


rule init_DESeq2:
    """ Create a RDS file for the different variables """
    input:
        counts = rules.prepare_data_for_DESeq2.output.counts,
        samples = rules.prepare_data_for_DESeq2.output.samples
    output:
        rds = "results/DESeq2/{ICAmodel}/deseq2_{variable}_base_{tool}.rds"
    log:
        "logs/DESeq2/init_DESeq2/{ICAmodel}_{variable}_{tool}.log"
    conda:
        "../envs/DESeq2.yaml"
    threads:
        40
    script:
        "../scripts/DEGs/2_init_DESeq2.R"


rule DESeq2:
    """ Calculate the DEA """
    input:
        rds = rules.init_DESeq2.output.rds
    output:
        results = "results/DESeq2/{ICAmodel}/{variable}/{tool}_vs_{tool2}.csv"
    conda:
        "../envs/DESeq2.yaml"
    log:
        "logs/DESeq2/DESeq2_{ICAmodel}/{variable}_{tool}_{tool2}.log"
    threads:
        40
    script:
        "../scripts/DEGs/3_DESeq2.R"


rule DESeq2_volcano_plot:
    """
        Volcano plots of the different workflow variables, where genes above a
        certain fold_change and p-value tresholds are colored.
    """
    input:
        DEGs = lambda w: get_DEG_results("{ICAmodel}".format(**w), config)
    output:
        plot = "results/DESeq2/plots/{ICAmodel}.png"
    params:
        fold_change = 2,
        pvalue = 1e-35
    conda:
        "../envs/ICA_python.yaml"
    script:
        "../scripts/DEGs/4_DESeq2_volcano_plot.py"


rule DESeq2_volcano_plot_comp:
    """
        Volcano plots of the different workflow variables, where significant
        negative and positives genes are colored.
    """
    input:
        DEGs = lambda w: get_DEG_results("{ICAmodel}".format(**w), config),
        gene_list = "results/ICA/sklearnFastICA/{ICA_path}/sigma_{sigma}/gene_list/comp_{comp}.txt",
        HGNC = rna_seq_cartesian_product("data/references/hgnc.txt")
    output:
        plot = "results/DESeq2/plots/{ICA_path}/{ICAmodel}_sigma{sigma}_comp{comp}.svg"
    conda:
        "../envs/ICA_python.yaml"
    script:
        "../scripts/DEGs/5_DESeq2_volcano_plot_comp.py"
