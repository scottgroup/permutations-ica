configfile: "config.json"

# Defining wildcards constraints
wildcard_constraints:
    no_constraint = ""


# Including rules
include: "rules/running_ICA.smk"


rule all:
    input:
        expand(
            "results/ICA/sklearnFastICA/counts_{isNaN}/M{M}_n10/dendrogram.png",
            M=range(7,25),
            isNaN=['noNaN']
        )
