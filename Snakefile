configfile: "config.json"

# Defining wildcards constraints
wildcard_constraints:
    no_constraint = ""


# Including rules
include: "rules/running_pyProDenICA.smk"


rule all:
    input:
        "results/ICA/pyProDenICA/counts_toy/M6_n10/components.tsv"
