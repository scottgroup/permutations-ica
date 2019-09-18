
rule running_ICA:
    input:
        dataset = "data/counts_noNaN.tsv"
    output:
        raw_components = "results/ICA/{ICAmethod}/{dataset}/M{M}_n{n}_std{std}/raw_components.tsv",
        fit_min = "results/ICA/{ICAmethod}/{dataset}/M{M}_n{n}_std{std}/fit_min.txt"
    params:
        max_it = 50000,
        tolerance = 1e-18
    threads:
        32
    conda:
        "../envs/ICA_python.yaml"
    script:
        "../scripts/running_ICA/running_ICA.py"


rule flipping_ICA_components:
    input:
        components = rules.running_ICA.output.raw_components
    output:
        components = "results/ICA/{ICAmethod}/{dataset}/M{M}_n{n}_std{std}/components.tsv",
        corr_components =  "results/ICA/{ICAmethod}/{dataset}/M{M}_n{n}_std{std}/corr_components.json",
    params:
        minimal_occurence = 0.20,
        threshold = 0.90
    conda:
        "../envs/ICA_python.yaml"
    script:
        "../scripts/running_ICA/flipping_ICA_components.py"


rule ICA_components_dendrogram:
    input:
        components = rules.flipping_ICA_components.output.components
    output:
        plot = "results/ICA/{ICAmethod}/{dataset}/M{M}_n{n}_std{std}/dendrogram.png",
        plot2 = "results/ICA/{ICAmethod}/{dataset}/plots/M{M}_n{n}_std{std}.png"
    conda:
        "../envs/ICA_python.yaml"
    script:
        "../scripts/plotting_ICA/plotting_dendrogram.py"
