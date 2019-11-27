
def get_components(wildcards):
    return expand(
        "results/ICA/{{ICAmethod}}/{{dataset}}/M{M}_n{{n}}_std{{std}}/components.tsv",
        M=range(int(wildcards.min), int(wildcards.max)+1)
    )



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
        plot = directory("results/ICA/{ICAmethod}/{dataset}/{ICA_run}/sigma_{sigma}/projection")
    params:
        fpath = lambda wildcards: "results/ICA/{ICAmethod}/{dataset}/{ICA_run}/sigma_{sigma}/projection/comp_{{comp}}.svg"
    conda:
        "../envs/ICA_python.yaml"
    script:
        "../scripts/plotting_ICA/plotting_component_projections.py"


rule plotting_M_stability_data:
    """
        Getting data for M_stability
    """
    input:
        get_components
    output:
        data = "results/ICA/{ICAmethod}/{dataset}/combine_{min}to{max}_n{n}_std{std}/M_stability.tsv"
    conda:
        "../envs/ICA_python.yaml"
    script:
        "../scripts/plotting_ICA/plotting_M_stability_data.py"


rule plotting_M_stability:
    """
        Plotting block_pearson score for every M
    """
    input:
        data = rules.plotting_M_stability_data.output.data
    output:
        plot = "results/ICA/{ICAmethod}/{dataset}/combine_{min}to{max}_n{n}_std{std}/M_stability.svg"
    conda:
        "../envs/ICA_python.yaml"
    script:
        "../scripts/plotting_ICA/plotting_M_stability.py"



rule plotting_distribution_grid:
    """
    """
    input:
        proj = "results/ICA/{ICAmethod}/{dataset}/{ICA_run}/filtered_components/sigma_{sigma}/projection.tsv",
        comp_list = "results/ICA/{ICAmethod}/{dataset}/{ICA_run}/filtered_components/sigma_{sigma}/comp_list.txt"
    output:
        "results/ICA/{ICAmethod}/{dataset}/{ICA_run}/sigma_{sigma}/grid.svg"
    conda:
        "../envs/ICA_python.yaml"
    script:
        "../scripts/plotting_ICA/plotting_distribution_grid.py"


rule plotting_heatmap_components:
    """
    """
    input:
        proj = "results/ICA/{ICAmethod}/{dataset}/{ICA_run}/filtered_components/sigma_{sigma}/projection.tsv",
        comp_list = "results/ICA/{ICAmethod}/{dataset}/{ICA_run}/filtered_components/sigma_{sigma}/comp_list.txt"
    output:
        plot = "results/ICA/{ICAmethod}/{dataset}/{ICA_run}/sigma_{sigma}/heatmap_components.svg"
    params:
        k_neighbour = 50
    conda:
        "../envs/ICA_python.yaml"
    script:
        "../scripts/plotting_ICA/plotting_heatmap_components.py"
