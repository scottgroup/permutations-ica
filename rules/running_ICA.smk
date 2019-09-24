

def get_components_range(wildcards):
    return expand(
        "results/ICA/{{ICAmethod}}/{{dataset}}/M{M}_n{{n}}_std{{std}}/components_mean.tsv",
        M=range(int(wildcards.min), int(wildcards.max)+1)
    )



rule running_ICA:
    """
        Running an ICA model on the data

        ICAmethod -> Specifies which model is being used. All methods are
            hosted in the scripts/running_ICA/ICAmethods folder and they all
            share the same Dataset object and preprocessing.
        dataset -> Specifies the subset of data chosen. All datasets must be
            defined in config.json under ICA_datasets as a dictionary where
            each key is a variable, and each value a list of possible variable
            states.
        M -> Number of components to generate.
        n -> Number of iteration. Each iteration only differs by the
            optimisation starting point.
    """
    input:
        counts = lambda w: config['ICA_datasets']['{dataset}'.format(**w)]['counts']
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
        "../scripts/running_ICA/1_running_ICA.py"


rule flipping_ICA_components:
    """
        Takes raw output from rules.running_ICA.running_ICA

        Since each ICA run is independent, correlation are calculated to find
        similar components from each run. Components are also "flipped" if they
        have a strong negative correlation.
    """
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
        "../scripts/running_ICA/2_flipping_ICA_components.py"


rule ICA_components_dendrogram:
    """
        Plots the correlation between all components from all iterations of an
        ICA model.
    """
    input:
        components = rules.flipping_ICA_components.output.components
    output:
        plot = "results/ICA/{ICAmethod}/{dataset}/M{M}_n{n}_std{std}/dendrogram.png",
        plot2 = "results/ICA/{ICAmethod}/{dataset}/plots/M{M}_n{n}_std{std}.png"
    conda:
        "../envs/ICA_python.yaml"
    script:
        "../scripts/plotting_ICA/plotting_dendrogram.py"


rule merging_n_components:
    """
        Calculating mean and std of the n correlated components to merge them
        into single components.
    """
    input:
        components = "results/{ICA_path}/components.tsv",
        corr_components = "results/{ICA_path}/corr_components.json"
    output:
        components_mean = "results/{ICA_path}/components_mean.tsv",
        components_std = "results/{ICA_path}/components_std.tsv"
    conda:
        "../envs/ICA_python.yaml"
    script:
        "../scripts/running_ICA/3_merging_n_components.py"


rule filter_sigma_components:
    """
        Assuming a normal distribution of the genes' weight for each
        components, keeping only gene that are sigma sigma away from mean.
    """
    input:
        components_mean = rules.merging_n_components.output.components_mean
    output:
        filt_genes = "results/{ICA_path}/filtered_components/sigma_{sigma}/components.tsv"
    conda:
        "../envs/ICA_python.yaml"
    script:
        "../scripts/running_ICA/4_filter_sigma_components.py"


rule dataset_projection_on_filtered_components:
    """

    """
    input:
        counts = lambda w: config['ICA_datasets']['{dataset}'.format(**w)]['counts'],
        components = "results/ICA/{ICAmethod}/{dataset}/{ICA_run}/filtered_components/sigma_{sigma}/components.tsv"
    output:
        projection = "results/ICA/{ICAmethod}/{dataset}/{ICA_run}/filtered_components/sigma_{sigma}/projection.tsv"
    conda:
        "../envs/ICA_python.yaml"
    script:
        "../scripts/running_ICA/5_dataset_projection_on_filt_comps.py"


rule component_agnostic_model:
    """

    """
    input:
        get_components_range
    output:
        components = "results/ICA/{ICAmethod}/{dataset}/combine_{min}to{max}_n{n}_std{std}/components.tsv",
        corr_components = "results/ICA/{ICAmethod}/{dataset}/combine_{min}to{max}_n{n}_std{std}/correlated_components.json"
    params:
        threshold = 0.95
    conda:
        "../envs/ICA_python.yaml"
    script:
        "../scripts/running_ICA/6_component_agnostic_model.py"

