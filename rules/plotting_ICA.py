from snakemake.io import expand

def get_ICA_plotting(configs, datasets):
    """
    Return all the different plots linked to the ICA models.
    """
    files = list()

    for dataset in datasets:

        # Loading config
        config = configs['ICA_datasets'][dataset]['params']
        config['dataset'] = dataset

        # Adding M list to params
        config['M'] = list(range(config['min'], config['max']+1))

        # Calculating all ICA_run
        ICA_run = list()
        ICA_run.extend(expand("M{M}_n{n}_std{std}", **config))
        ICA_run.extend(expand("combine_{min}to{max}_n{n}_std{std}", **config))
        config['ICA_run'] = ICA_run

        # Adding dendrogram and correlation matrix for the components
        str = "results/ICA/{ICAmethod}/{dataset}/{ICA_run}/{plot}.png"
        files.extend(expand(str, plot=['dendrogram', 'corr'], **config))

        # Adding component projections on variables
        str = "results/ICA/{ICAmethod}/{dataset}/{ICA_run}/sigma_{sigma}/projection"
        files.extend(expand(str, **config))

        # Adding M_stability plot
        str = "results/ICA/{ICAmethod}/{dataset}/combine_{min}to{max}_n{n}_std{std}/M_stability.svg"
        files.extend(expand(str, **config))

        # Adding heatmaps
        str = "results/ICA/{ICAmethod}/{dataset}/{ICA_run}/sigma_{sigma}/heatmap_components.svg"
        files.extend(expand(str, **config))

    return files
