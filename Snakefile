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


# Adding subworkflow
include: "subworkflows/subworkflows.smk"


# Including rules
include: "rules/running_ICA.smk"
include: "rules/plotting_ICA.smk"
include: "rules/analyse_ICA.smk"
include: "rules/DEGs.smk"
include: "rules/pseudogene_analysis.smk"


# Defining model to run
ICAruns = list(config['ICA_datasets'].keys())


rule all:
    input:
        # pseudogene_parent("results/blat_score.tsv"),
        rna_seq_cartesian_product("results/cartesian_product/tissues_NaN.tsv"),
        # rna_seq_cartesian_product("results/rnaseq/geneCoverage/EEF1A1_plot.png"),
        # running_ICA.get_ICA_running(config, ICAruns),
        # plotting_ICA.get_ICA_plotting(config, ICAruns),
        # DEGs.get_DESeq(config, ICAruns)
