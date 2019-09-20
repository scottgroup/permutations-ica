configfile: "config.json"

# Defining wildcards constraints
wildcard_constraints:
    dataset = "({})".format("|".join(config["ICA_datasets"].keys())),
    ICAmethod = "sklearnFastICA",
    std="[0-9]+"


# Including rules
include: "rules/running_ICA.smk"

# Defining variables
M = range(9, 16)
std = [0, 1, 2, 3]
datasets = ["counts_noNaN", "counts_NaN"]
sigma = [1, 4, 9]


rule all:
    input:
        expand(
            "results/ICA/sklearnFastICA/{dataset}/M{M}_n10_std{std}/filtered_components/sigma_{sigma}/projection.tsv",
            M=M, std=std, dataset=datasets, sigma=sigma
        ),
        expand(
            "results/ICA/sklearnFastICA/{dataset}/M{M}_n10_std{std}/dendrogram.png",
            M=M, std=std, dataset=datasets, sigma=sigma
        )
