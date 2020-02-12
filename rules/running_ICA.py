from snakemake.io import expand


def get_components_range(wildcards):
    """ Returns list of paths with wildards. M is specified, from min to max """
    return expand(
        "results/ICA/{{ICAmethod}}/{{ICAmodel}}/M{M}_n{{n}}_std{{std}}/components_mean.tsv",
        M=range(int(wildcards.min), int(wildcards.max)+1)
    )


def get_ICA_running(configs, ICAmodels):
    """
    For specified M in ICAmodel parameters, run the different analysis scripts
    """
    files = list()

    for ICAmodel in ICAmodels:

        # Loading config
        config = configs['ICA_models'][ICAmodel]['params']
        config['ICAmodel'] = ICAmodel

        # Adding M list to params
        config['M'] = config['to_analyse']

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
        path = "results/ICA/{ICAmethod}/{ICAmodel}/{ICA_run}/sigma_{sigma}/filtered_components/projection.tsv"
        files.extend(expand(path, **config))

    return files
