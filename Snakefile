configfile: "config.json"

import rules.running_ICA as running_ICA
import rules.plotting_ICA as plotting_ICA
import rules.DEGs as DEGs

# Defining wildcards constraints
wildcard_constraints:
    dataset = "({})".format("|".join(config["ICA_datasets"].keys())),
    ICAmethod = "({})".format("|".join(["sklearnFastICA", "consICA"])),
    std="[0-9]+",
    sigma="[0-9]+"

# Including rules
include: "rules/running_ICA.smk"
include: "rules/plotting_ICA.smk"
include: "rules/analyse_ICA.smk"
include: "rules/DEGs.smk"
include: "rules/pseudogene_analysis.smk"

# Defining model to run
ICAruns = list(config['ICA_datasets'].keys())

# Adding subworkflow
subworkflow pseudogene_parent:
    workdir:
        "subworkflows/pseudogene_parent"
    snakefile:
        "subworkflows/pseudogene_parent/Snakefile"
    configfile:
        "subworkflows/pseudogene_parent/config.json"


rule all:
    input:
        # pseudogene_parent("results/blat_score.tsv")
        running_ICA.get_ICA_running(config, ICAruns),
        plotting_ICA.get_ICA_plotting(config, ICAruns),
        # DEGs.get_DESeq(config, ICAruns)
