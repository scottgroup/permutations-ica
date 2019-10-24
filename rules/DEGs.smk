from DEGs import get_DEG_inputs


rule prepare_data_for_DESeq2:
    input:
        dataset = lambda w: config['ICA_datasets']['{dataset}'.format(**w)]['params']['counts'],
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
        32
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
        32
    script:
        "../scripts/DEGs/3_DESeq2.R"


rule analyze_DESeq2:
    """ TODO : rework this function """
    input:
        DEGs = lambda w: get_DEG_inputs("{dataset}".format(**w), config),
        # hgnc = "data/references/hgnc.txt"
    output:
        "something_{dataset}"
    params:
        fold_change = 1,
        pvalue = 0.00000001
    conda:
        "../envs/ICA_python.yaml"
    script:
        "../scripts/DEGs/analyze_DESeq2.py"
