from snakemake.io import expand

def get_ICA_running(config, dataset):
    """
    Return all necessary files for computations of the model.

    This function intentionally exports too many files. Some of them are not
    necessary for a sound completion of the model. The redudant file names are
    intented to give modularity.
    """
    files = list()

    # Loading config
    config = config['ICA_datasets'][dataset]['params']
    config['dataset'] = dataset

    # Adding M list to params
    config['M'] = list(range(config['min'], config['max']+1))

    # Calculating all ICA_run
    ICA_run = list()
    ICA_run.extend(expand("M{M}_n{n}_std{std}", **config))
    ICA_run.extend(expand("combine_{min}to{max}_n{n}_std{std}", **config))
    config['ICA_run'] = ICA_run

    # ICA model for a specific M
    path = "results/ICA/{ICAmethod}/{dataset}/M{M}_n{n}_std{std}/components.tsv"
    files.extend(expand(path, **config))

    # ICA model for combined M
    path = "results/ICA/{ICAmethod}/{dataset}/combine_{min}to{max}_n{n}_std{std}/components.tsv"
    files.extend(expand(path, **config))

    # ICA projections of components
    path = "results/ICA/{ICAmethod}/{dataset}/{ICA_run}/filtered_components/sigma_{sigma}/projection.tsv"
    files.extend(expand(path, **config))

    return files
