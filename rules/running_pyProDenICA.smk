
rule running_ICA:
    input:
        dataset = "data/counts.tsv"
    output:
        components = "results/ICA/pyProDenICA/{dataset}/M{M}_n{n}/components.tsv",
        fit_min = "results/ICA/pyProDenICA/{dataset}/M{M}_n{n}/fit_min.txt"
    params:
        max_it = 1000
    conda:
        "../envs/ICA_python.yaml"
    script:
        "../scripts/running_pyProDenICA/running_pyProDenICA.py"
