

rule ICA_components_dendrogram:
    """
        Plots the correlation between all components from all iterations of an
        ICA model, using a dendrogram.
        TODO: Filter when there is no acceptable component
    """
    input:
        components = "results/ICA/{ICAmethod}/{dataset}/{ICA_run}/components.tsv"
    output:
        plot = "results/ICA/{ICAmethod}/{dataset}/{ICA_run}/dendrogram.png",
        plot2 = "results/ICA/{ICAmethod}/{dataset}/plots/{ICA_run}_dendrogram.png"
    params:
        threshold = 0.90
    conda:
        "../envs/ICA_python.yaml"
    script:
        "../scripts/plotting_ICA/plotting_dendrogram.py"


rule ICA_components_corr:
    """
        Plots the correlation between all components from all iterations of an
        ICA model, using a correlation matrix.
        TODO: Filter when there is no acceptable component
    """
    input:
        components = "results/ICA/{ICAmethod}/{dataset}/{ICA_run}/components.tsv"
    output:
        plot = "results/ICA/{ICAmethod}/{dataset}/{ICA_run}/corr.png",
        plot2 = "results/ICA/{ICAmethod}/{dataset}/plots/{ICA_run}_corr.png"
    conda:
        "../envs/ICA_python.yaml"
    script:
        "../scripts/plotting_ICA/ICA_components_corr.py"
