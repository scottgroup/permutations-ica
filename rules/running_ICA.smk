

def get_components_range(wildcards):
    return expand(
        "results/ICA/{{ICAmethod}}/{{dataset}}/M{M}_n{{n}}_std{{std}}/components_mean.tsv",
        M=range(int(wildcards.min), int(wildcards.max)+1)
    )


def get_counts(wildcards):
    """ """
    fname = config['ICA_datasets'][wildcards.dataset]['params']['counts']
    return rna_seq_cartesian_product("results/cartesian_product/{dataset}.tsv".format(dataset=fname)) 


rule running_sklearnFastICA:
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
        counts = get_counts,
    output:
        components = temp("results/ICA/sklearnFastICA/{dataset}/M{M}_n{n}_std{std}/raw_components.tsv"),
        fit_min = "results/ICA/sklearnFastICA/{dataset}/M{M}_n{n}_std{std}/fit_min.txt",
    params:
        max_it = 50000,
        tolerance = 1e-20
    threads:
        40
    conda:
        "../envs/ICA_python.yaml"
    script:
        "../scripts/running_ICA/1_running_sklearnFastICA.py"


rule flipping_ICA_components:
    """
        Takes raw output from rules.running_ICA.running_ICA

        Since each ICA run is independent, correlation are calculated to find
        similar components from each run. Components are also "flipped" if they
        have a strong negative correlation.
    """
    input:
        components = "results/ICA/{ICAmethod}/{dataset}/M{M}_n{n}_std{std}/raw_components.tsv",
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


rule merging_n_components:
    """
        Calculating mean and std of the n correlated components to merge them
        into single components.
    """
    input:
        components = "results/ICA/{ICAmethod}/{dataset}/{ICA_run}/components.tsv",
        corr_components = "results/ICA/{ICAmethod}/{dataset}/{ICA_run}/corr_components.json",
    output:
        components_mean = "results/ICA/{ICAmethod}/{dataset}/{ICA_run}/components_mean.tsv",
        components_std = "results/ICA/{ICAmethod}/{dataset}/{ICA_run}/components_std.tsv",
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
        components_mean = "results/ICA/{ICAmethod}/{dataset}/{ICA_run}/components_mean.tsv",
    output:
        filt_genes = "results/ICA/{ICAmethod}/{dataset}/{ICA_run}/filtered_components/sigma_{sigma}/components.tsv",
    conda:
        "../envs/ICA_python.yaml"
    script:
        "../scripts/running_ICA/4_filter_sigma_components.py"


rule component_agnostic_model:
    """

    """
    input:
        get_components_range
    output:
        components = "results/ICA/{ICAmethod}/{dataset}/combine_{min}to{max}_n{n}_std{std}/components.tsv",
        corr_components = "results/ICA/{ICAmethod}/{dataset}/combine_{min}to{max}_n{n}_std{std}/corr_components.json",
    params:
        threshold = 0.90
    conda:
        "../envs/ICA_python.yaml"
    script:
        "../scripts/running_ICA/5_component_agnostic_model.py"


rule dataset_projection_on_filtered_components:
    """

    """
    input:
        counts = get_counts,
        components = "results/ICA/{ICAmethod}/{dataset}/{ICA_run}/filtered_components/sigma_{sigma}/components.tsv",
    output:
        projection = "results/ICA/{ICAmethod}/{dataset}/{ICA_run}/filtered_components/sigma_{sigma}/projection.tsv",
    conda:
        "../envs/ICA_python.yaml"
    script:
        "../scripts/running_ICA/6_dataset_projection_on_filt_comps.py"


rule dataset_variable_boolean:
    input:
        counts = get_counts,
    output:
        var_bool = "results/ICA/variable_boolean/{dataset}.tsv",
    conda:
        "../envs/ICA_python.yaml"
    script:
        "../scripts/running_ICA/7_dataset_variable_boolean.py"


rule projection_boolean_correlation:
    input:
        projection = rules.dataset_projection_on_filtered_components.output.projection,
        var_bool = rules.dataset_variable_boolean.output.var_bool,
    output:
        proj_corr = "results/ICA/{ICAmethod}/{dataset}/{ICA_run}/filtered_components/sigma_{sigma}/projection_correlation.tsv",
    conda:
        "../envs/ICA_python.yaml"
    script:
        "../scripts/running_ICA/8_projection_boolean_correlation.py"


rule list_distribution_component:
    """
        TODO : generates NaN on pearson. Okay? Need to filter?
    """
    input:
        proj_corr = rules.projection_boolean_correlation.output.proj_corr,
        var_bool = rules.dataset_variable_boolean.output.var_bool,
    output:
        list_comp = "results/ICA/{ICAmethod}/{dataset}/{ICA_run}/filtered_components/sigma_{sigma}/comp_list.txt",
    conda:
        "../envs/ICA_python.yaml"
    script:
        "../scripts/running_ICA/9_list_distribution_component.py"


rule extracting_gene_list:
    input:
        filt_genes = rules.filter_sigma_components.output.filt_genes
    output:
        tkn = "results/ICA/{ICAmethod}/{dataset}/{ICA_run}/gene_list/comp_sigma{sigma}/.tkn"
    params:
        directory = "results/ICA/{ICAmethod}/{dataset}/{ICA_run}/gene_list/comp_sigma{sigma}"
    conda:
        "../envs/ICA_python.yaml"
    script:
        "../scripts/running_ICA/10_extracting_gene_list.py"


rule running_sklearnFastICA_bootstraped:
    """
        Same as rule.running_sklearnFastICA.

        Except : Adding a bootstrapped parameter, generating components using
        {boot}% of original data.
    """
    input:
        counts = get_counts,
    output:
        components = "results/ICA/sklearnFastICA/{dataset}/M{M}_n{n}_std{std}/stability/components_{boot}_strapped.tsv",
        fit_min = "results/ICA/sklearnFastICA/{dataset}/M{M}_n{n}_std{std}/stability/fit_min_{boot}_strapped.txt",
    params:
        max_it = 50000,
        tolerance = 1e-20
    threads:
        32
    conda:
        "../envs/ICA_python.yaml"
    script:
        "../scripts/running_ICA/1_running_sklearnFastICA.py"


rule bootstrapped_stability_correlation:
    """

    """
    input:
        components = "results/ICA/{ICAmethod}/{dataset}/M{M}_n{n}_std{std}/components_mean.tsv",
        boot_components = "results/ICA/{ICAmethod}/{dataset}/M{M}_n{n}_std{std}/stability/components_{boot}_strapped.tsv",
    output:
        plot = "results/ICA/{ICAmethod}/{dataset}/M{M}_n{n}_std{std}/stability/{boot}_analysis.png",
    params:
        min_correlation = 0.75
    conda:
        "../envs/ICA_python.yaml"
    script:
        "../scripts/running_ICA/12_bootstrapped_stability_correlation.py"


rule bootstrapped_stability_n_genes:
    """

    """
    input:
        components = "results/ICA/{ICAmethod}/{dataset}/M{M}_n{n}_std{std}/components_mean.tsv",
        boot_components = "results/ICA/{ICAmethod}/{dataset}/M{M}_n{n}_std{std}/stability/components_{boot}_strapped.tsv",
    output:
        plot = "results/ICA/{ICAmethod}/{dataset}/M{M}_n{n}_std{std}/stability/{boot}_analysis_n_genes.png",
    params:
        sigma = 4,
        min_genes = 0.50
    conda:
        "../envs/ICA_python.yaml"
    script:
        "../scripts/running_ICA/12_bootstrapped_stability_n_genes.py"
