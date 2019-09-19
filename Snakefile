configfile: "config.json"

# Defining wildcards constraints
wildcard_constraints:
    dataset = "({})".format("|".join(config["ICA_datasets"].keys())),
    ICAmethod = "sklearnFastICA",
    std="[0-9]+"


# Including rules
include: "rules/running_ICA.smk"


rule all:
    input:
        # expand(
        #     "results/ICA/sklearnFastICA/counts_{isNaN}/plots/M{M}_n20_std{std}.png",
	    #     M=range(9,16),
        #     isNaN=['noNaN', 'NaN'],
	    #     std=[0, 1, 2, 3, 4]
        # )
        expand(
            "results/ICA/sklearnFastICA/{dataset}/M{M}_n10_std{std}/components_mean.tsv",
            M=range(4,6),
            std=[3],
            dataset=["counts_noNaN"]
        ),
        expand(
            "results/ICA/sklearnFastICA/{dataset}/M{M}_n10_std{std}/filtered_components/sigma_{sigma}.tsv",
            M=range(4,6),
            std=[3],
            dataset=["counts_noNaN"],
            sigma=[1, 4, 9]
        ),
