configfile: "config.json"

# Defining wildcards constraints
wildcard_constraints:
    dataset = "({})".format("|".join(config["ICA_datasets"].keys())),
    ICAmethod = "sklearnFastICA"


# Including rules
include: "rules/running_ICA.smk"


rule all:
    input:
        expand(
            # "results/ICA/sklearnFastICA/counts_{isNaN}/M{M}_n20_std{std}/dendrogram.png",
            "results/ICA/sklearnFastICA/counts_{isNaN}/plots/M{M}_n20_std{std}.png",
	    M=range(9,16),
            isNaN=['noNaN', 'NaN'],
	    std=[0, 1, 2, 3, 4]
        )
