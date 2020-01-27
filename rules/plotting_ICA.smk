# Paths
ICAmodel_path = config["path"]["ICAmodel"]
combine_path = config["path"]["combineModel"]


def get_components(wildcards):
    return expand(
        "results/ICA/{{ICAmethod}}/{{ICAmodel}}/M{M}_n{{n}}_std{{std}}/components.tsv",
        M=range(int(wildcards.min), int(wildcards.max)+1)
    )


rule plotting_components_dendrogram:
    """
        Plots the correlation between all components from all iterations of an
        ICA model, using a dendrogram.
    """
    input:
        components = ICAmodel_path + "components.tsv"
    output:
        plot = ICAmodel_path + "dendrogram.svg",
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
    """
    input:
        components = ICAmodel_path + "components.tsv"
    output:
        plot = ICAmodel_path + "corr.svg",
    conda:
        "../envs/ICA_python.yaml"
    script:
        "../scripts/plotting_ICA/plotting_components_corr.py"


rule plotting_component_projections:
    """
        Plots density graphs for the different pipeline variables, with each
        option colored differently.

        Output is a directory, creating projections for all components at once.
    """
    input:
        projection = rules.dataset_projection_on_filtered_components.output.projection
    output:
        plot = directory(ICAmodel_path + "sigma_{sigma}/projection")
    params:
        fpath = lambda wildcards: ICAmodel_path + "sigma_{sigma}/projection/comp_{{comp}}.svg"
    conda:
        "../envs/ICA_python.yaml"
    script:
        "../scripts/plotting_ICA/plotting_component_projections.py"


rule getting_M_stability_data:
    """
        Getting data for M_stability
    """
    input:
        get_components
    output:
        data = combine_path + "M_stability.tsv"
    conda:
        "../envs/ICA_python.yaml"
    script:
        "../scripts/plotting_ICA/getting_M_stability_data.py"


rule plotting_M_stability:
    """
        Plotting block_pearson score for every M
    """
    input:
        data = rules.getting_M_stability_data.output.data
    output:
        plot = combine_path + "M_stability.svg"
    conda:
        "../envs/ICA_python.yaml"
    script:
        "../scripts/plotting_ICA/plotting_M_stability.py"



rule plotting_distribution_grid:
    """
        Pairwise distributions of selected components. The components must be
        indicated as a list in the params. The components number are 0-based in
        the input and 1-based in the output.
    """
    input:
        proj = ICAmodel_path + "filtered_components/sigma_{sigma}/projection.tsv",
        comp_list = ICAmodel_path + "filtered_components/sigma_{sigma}/comp_list.txt"
    output:
        plot = ICAmodel_path + "sigma_{sigma}/grid.svg"
    params:
        comps = ['8',  '9', '12', '14']
    conda:
        "../envs/ICA_python.yaml"
    script:
        "../scripts/plotting_ICA/plotting_distribution_grid.py"


rule plotting_heatmap_components:
    """
        Create the KNN score heatmap for correlation between expression modes
        and methodological variables.
    """
    input:
        proj = ICAmodel_path + "filtered_components/sigma_{sigma}/projection.tsv",
        comp_list = ICAmodel_path + "filtered_components/sigma_{sigma}/comp_list.txt"
    output:
        plot = ICAmodel_path + "sigma_{sigma}/heatmap_components.svg"
    params:
        k_neighbour = 50
    conda:
        "../envs/ICA_python.yaml"
    script:
        "../scripts/plotting_ICA/plotting_heatmap_components.py"
