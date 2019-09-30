
def get_plots():
    pass


rule plotting_components_dendrogram:
    """
        Plots the correlation between all components from all iterations of an
        ICA model, using a dendrogram.
        TODO: Filter when there is no acceptable component
    """
    input:
        components = "results/ICA/{ICAmethod}/{dataset}/{ICA_run}/components.tsv"
    output:
        plot = "results/ICA/{ICAmethod}/{dataset}/{ICA_run}/dendrogram.png",
        plot2 = "results/ICA/{ICAmethod}/{dataset}/plots/{ICA_run}/dendrogram.png"
    params:
        threshold = 0.90
    conda:
        "../envs/ICA_python.yaml"
    script:
        "../scripts/plotting_ICA/plotting_components_dendrogram.py"


rule plotting_components_corr:
    """
        Plots the correlation between all components from all iterations of an
        ICA model, using a correlation matrix.
        TODO: Filter when there is no acceptable component
    """
    input:
        components = "results/ICA/{ICAmethod}/{dataset}/{ICA_run}/components.tsv"
    output:
        plot = "results/ICA/{ICAmethod}/{dataset}/{ICA_run}/corr.png",
        plot2 = "results/ICA/{ICAmethod}/{dataset}/plots/{ICA_run}/corr.png"
    conda:
        "../envs/ICA_python.yaml"
    script:
        "../scripts/plotting_ICA/plotting_components_corr.py"


rule plotting_component_projections:
    """

    """
    input:
        projection = rules.dataset_projection_on_filtered_components.output.projection
    output:
        plot = directory("results/ICA/{ICAmethod}/{dataset}/plots/{ICA_run}/sigma_{sigma}/projection")
    params:
        fpath = lambda wildcards: "results/ICA/{ICAmethod}/{dataset}/plots/{ICA_run}/sigma_{sigma}/projection/comp_{{comp}}.png"
    conda:
        "../envs/ICA_python.yaml"
    script:
        "../scripts/plotting_ICA/plotting_component_projections.py"
