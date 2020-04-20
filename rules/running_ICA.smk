from running_ICA import get_components_range

ICAmodel_path = config["path"]["ICAmodel"]

def get_counts(wildcards):
    """ Returns path of the count file """
    fname = config['ICA_models'][wildcards.ICAmodel]['params']['counts']
    dataset = "results/cartesian_product/{dataset}.tsv".format(dataset=fname)
    if 'stranded' in wildcards.ICAmodel:
        return rna_seq_cartesian_product_replicate(dataset)
    else:
        return rna_seq_cartesian_product(dataset)


rule running_sklearnFastICA:
    """
        Running an ICA model on the data

        ICAmethod -> Specifies which model is being used. All methods are
            hosted in the scripts/running_ICA/ICAmethods folder and they all
            share the same Dataset object and preprocessing.
        dataset -> Specifies the subset of data chosen. All datasets must be
            defined in config.json under ICA_models as a dictionary where
            each key is a variable, and each value a list of possible variable
            states.
        M -> Number of components to generate.
        n -> Number of iteration. Each iteration only differs by the
            optimisation starting point.
        std -> Building the ICA model only with genes that have a standard
            deviation (of expression value) this is at least std*(mean standard
            deviation for all genes). Keeps all the gene at std = 0.
    """
    input:
        counts = get_counts
    output:
        components = temp(ICAmodel_path + "raw_components.tsv"),
        fit_min = ICAmodel_path + "fit_min.txt",
    params:
        max_it = 100000,
        tolerance = 1e-18
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
        components = rules.running_sklearnFastICA.output.components
    output:
        components = ICAmodel_path + "components.tsv",
        corr_components =  ICAmodel_path + "corr_components.json",
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
        components = rules.flipping_ICA_components.output.components,
        corr_components = rules.flipping_ICA_components.output.corr_components,
    output:
        components_mean = ICAmodel_path + "components_mean.tsv",
        components_std = ICAmodel_path + "components_std.tsv",
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
        filt_genes = ICAmodel_path + "sigma_{sigma}/filtered_components/components.tsv",
    conda:
        "../envs/ICA_python.yaml"
    script:
        "../scripts/running_ICA/4_filter_sigma_components.py"


rule dataset_projection_on_filtered_components:
    """
        Calculate the dot product between gene weights and gene expression data
    """
    input:
        counts = get_counts,
        components = rules.filter_sigma_components.output.filt_genes
    output:
        projection = ICAmodel_path + "sigma_{sigma}/filtered_components/projection.tsv",
    conda:
        "../envs/ICA_python.yaml"
    script:
        "../scripts/running_ICA/6_dataset_projection_on_filt_comps.py"


rule dataset_variable_boolean:
    """
        Creates a binary matrix with information of software used for each
        pipeline
    """
    input:
        counts = get_counts
    output:
        var_bool = "results/ICA/variable_boolean/{ICAmodel}.tsv",
    conda:
        "../envs/ICA_python.yaml"
    script:
        "../scripts/running_ICA/7_dataset_variable_boolean.py"


rule projection_boolean_correlation:
    """

    """
    input:
        projection = rules.dataset_projection_on_filtered_components.output.projection,
        var_bool = rules.dataset_variable_boolean.output.var_bool,
    output:
        proj_corr = ICAmodel_path + "sigma_{sigma}/filtered_components/projection_correlation.tsv",
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
        list_comp = ICAmodel_path + "sigma_{sigma}/filtered_components/comp_list.txt",
    conda:
        "../envs/ICA_python.yaml"
    script:
        "../scripts/running_ICA/9_list_distribution_component.py"


rule extracting_gene_list:
    """
        For all components in a model, create a .txt file with Ensembl IDs of
        the significant genes (respectively to sigma), separated by positive
        and negative weight.

        This rule only specifies an output token in lieu of all the generated
        files. This token must be use in inputs of rules using gene lists.
    """
    input:
        filt_genes = rules.filter_sigma_components.output.filt_genes
    output:
        tkn = ICAmodel_path + "sigma_{sigma}/gene_list/.tkn"
    params:
        directory = ICAmodel_path + "sigma_{sigma}/gene_list"
    conda:
        "../envs/ICA_python.yaml"
    script:
        "../scripts/running_ICA/10_extracting_gene_list.py"
