
def get_permutations(k):
    for i in range(1, k):
        for j in range(i+1, k+1):
            yield i-1, j-1

def get_DEG_inputs(dataset, config):
    DEG_experiments = list()
    fpath = "results/DESeq2/{{dataset}}/{variable}/{tool}_vs_{tool2}.csv"

    # Getting variable dictionnary
    var_dict = config['tools']
    for k, v in config['ICA_datasets'][dataset]['variables'].items():
        var_dict[k] = v

    # Iterate through the variables
    for variable in var_dict.keys():
        tools = var_dict[variable]
        for i, j in get_permutations(len(tools)):
            _fpath = fpath.format(
                variable=variable, tool=tools[i], tool2=tools[j]
            )
            DEG_experiments.append(_fpath)

    # Iterate through the tissues
    tissues = list(config['datasets'].keys())
    for i, j in get_permutations(len(tissues)):
        _fpath = fpath.format(
            variable='tissue', tool=tissues[i], tool2=tissues[j]
        )
        DEG_experiments.append(_fpath)

    return DEG_experiments
