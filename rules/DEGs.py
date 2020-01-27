from snakemake.io import expand

def get_permutations(k):
    """ Get permutation for pairwise comparisons """
    for i in range(1, k):
        for j in range(i+1, k+1):
            yield i-1, j-1


def get_DEG_results(ICAmodel, config):
    """ For an ICAmodel, creates the different files """
    DEG_experiments = list()
    fpath = "results/DESeq2/{{ICAmodel}}/{variable}/{tool}_vs_{tool2}.csv"

    # Getting variable dictionnary
    var_dict = config['tools']
    for k, v in config['ICA_models'][ICAmodel]['variables'].items():
        var_dict[k] = v

    # Iterate through the variables
    for variable in var_dict.keys():
        tools = var_dict[variable]
        for i, j in get_permutations(len(tools)):
            _fpath = fpath.format(
                variable=variable, tool=tools[i], tool2=tools[j]
            )
            DEG_experiments.append(_fpath)

    return DEG_experiments


def get_DESeq(config, ICAmodels):
    """ Run the DESeq analyses for all ICA models """
    all_files = list()

    for ICAmodel in ICAmodels:
        DEGs_exps = [
            file.format(ICAmodel=ICAmodel)
            for file in get_DEG_results(ICAmodel, config)
        ]
        all_files.extend(DEGs_exps)

    return all_files
