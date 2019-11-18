from DEGs import get_DEG_results


def get_counts(wildcards):
    """ """
    fname = config['ICA_datasets'][wildcards.dataset]['params']['counts']
    return rna_seq_cartesian_product("results/cartesian_product/{dataset}.tsv".format(dataset=fname))


rule prepare_data_for_DESeq2:
    input:
        dataset = get_counts,
    output:
        counts = "results/DESeq2/{dataset}/counts.tsv",
        samples = "results/DESeq2/{dataset}/samples.tsv"
    conda:
        "../envs/ICA_python.yaml"
    script:
        "../scripts/DEGs/1_prepare_data_for_DESeq2.py"


rule init_DESeq2:
    input:
        counts = rules.prepare_data_for_DESeq2.output.counts,
        samples = rules.prepare_data_for_DESeq2.output.samples
    output:
        rds = "results/DESeq2/{dataset}/deseq2_{variable}_base_{tool}.rds"
    log:
        "logs/DESeq2/init_DESeq2/{dataset}_{variable}_{tool}.log"
    conda:
        "../envs/DESeq2.yaml"
    threads:
        40
    script:
        "../scripts/DEGs/2_init_DESeq2.R"


rule DESeq2:
    input:
        rds = rules.init_DESeq2.output.rds
    output:
        results = "results/DESeq2/{dataset}/{variable}/{tool}_vs_{tool2}.csv"
    conda:
        "../envs/DESeq2.yaml"
    log:
        "logs/DESeq2/DESeq2_{dataset}/{variable}_{tool}_{tool2}.log"
    threads:
        40
    script:
        "../scripts/DEGs/3_DESeq2.R"


rule DESeq2_volcano_plot:
    """ TODO : rework this function """
    input:
        DEGs = lambda w: get_DEG_results("{dataset}".format(**w), config)
    output:
        plot = "results/DESeq2/plots/{dataset}.png"
    params:
        fold_change = 2,
        pvalue = 1e-8
    conda:
        "../envs/ICA_python.yaml"
    script:
        "../scripts/DEGs/4_DESeq2_volcano_plot.py"
