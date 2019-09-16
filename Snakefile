configfile: "config.json"

# Defining wildcards constraints
wildcard_constraints:
    no_constraint = ""


# Including rules
include: "rules/running_ICA.smk"


rule all:
    input:
        # "results/ICA/sklearnFastICA/counts_toy/M9_n10/components.tsv",
        expand(
            "results/ICA/sklearnFastICA/counts_toy/M{M}_n10/dendrogram.png",
            M=range(4,18)
        )
