from snakemake.io import expand

def find_optimal_M(configs, ICAmodels):
    """
    Return plots needed to select optimal M
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
        config['ICA_run'] = ICA_run

        # Adding dendrogram and correlation matrix for the components
        str = "results/ICA/{ICAmethod}/{ICAmodel}/{ICA_run}/{plot}.svg"
        files.extend(expand(str, plot=['dendrogram', 'corr'], **config))

        # Adding M_stability plot
        str = "results/ICA/{ICAmethod}/{ICAmodel}/combine_{min}to{max}_n{n}_std{std}/M_stability.svg"
        files.extend(expand(str, **config))

    return files



def get_ICA_plotting(configs, ICAmodels):
    """
    For specified M in ICAmodel parameters, run the different analysis plots
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
        config['ICA_run'] = ICA_run

        # Adding dendrogram and correlation matrix for the components
        str = "results/ICA/{ICAmethod}/{ICAmodel}/{ICA_run}/{plot}.svg"
        files.extend(expand(str, plot=['dendrogram', 'corr'], **config))

        # Adding component projections on variables
        str = "results/ICA/{ICAmethod}/{ICAmodel}/{ICA_run}/sigma_{sigma}/projection"
        files.extend(expand(str, **config))

        # Adding M_stability plot
        str = "results/ICA/{ICAmethod}/{ICAmodel}/combine_{min}to{max}_n{n}_std{std}/M_stability.svg"
        files.extend(expand(str, **config))

        # Adding heatmaps
        str = "results/ICA/{ICAmethod}/{ICAmodel}/{ICA_run}/sigma_{sigma}/heatmap_components.svg"
        files.extend(expand(str, **config))

    return files
