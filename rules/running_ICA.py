from snakemake.io import expand


def get_components_range(wildcards):
    """ Returns list of paths with wildards. M is specified, from min to max """
    return expand(
        "results/ICA/{{ICAmethod}}/{{ICAmodel}}/M{M}_n{{n}}_std{{std}}/components_mean.tsv",
        M=range(int(wildcards.min), int(wildcards.max)+1)
    )


def get_ICA_running(configs, ICAmodels):
    """
    Return all necessary files for computations of the model.

    This function intentionally exports too many files. Some of them are not
    necessary for a sound completion of the model. The redudant file names are
    intented to give modularity.
    """
    files = list()

    for ICAmodel in ICAmodels:

        # Loading config
        config = configs['ICA_models'][ICAmodel]['params']
        config['ICAmodel'] = ICAmodel

        # Adding M list to params
        config['M'] = list(range(config['min'], config['max']+1))

        # Calculating all ICA_run
        ICA_run = list()
        ICA_run.extend(expand("M{M}_n{n}_std{std}", **config))
        ICA_run.extend(expand("combine_{min}to{max}_n{n}_std{std}", **config))
        config['ICA_run'] = ICA_run

        # ICA model for a specific M
        path = "results/ICA/{ICAmethod}/{ICAmodel}/M{M}_n{n}_std{std}/components.tsv"
        files.extend(expand(path, **config))

        # ICA model for combined M
        path = "results/ICA/{ICAmethod}/{ICAmodel}/combine_{min}to{max}_n{n}_std{std}/components.tsv"
        files.extend(expand(path, **config))

        # ICA projections of components
        path = "results/ICA/{ICAmethod}/{ICAmodel}/{ICA_run}/filtered_components/sigma_{sigma}/projection.tsv"
        files.extend(expand(path, **config))

    return files
