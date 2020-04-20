configfile: "config.json"

import rules.running_ICA as running_ICA
import rules.plotting_ICA as plotting_ICA
import rules.DEGs as DEGs


# Defining wildcards constraints
wildcard_constraints:
    dataset = "({})".format("|".join(config["ICA_models"].keys())),
    ICAmethod = "({})".format("|".join(["sklearnFastICA"])),
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
include: "rules/genomic_overlap_analysis.smk"

# Defining model to run
ICAmodels = list(config['ICA_models'].keys())


rule all:
    input:
        # Running the RNA-seq pipelines
        rna_seq_cartesian_product("results/cartesian_product/tissues_NaN.tsv"),
        rna_seq_cartesian_product("results/cartesian_product/tissues_noNaN.tsv"),
        # Running the alternative RNA-seq pipelines
        rna_seq_cartesian_product_replicate("results/cartesian_product/tissues_stranded_NaN.tsv"),
        rna_seq_cartesian_product_replicate("results/cartesian_product/tissues_unstranded_NaN.tsv"),
        # Generating the basic plots to find the optimal M
        plotting_ICA.find_optimal_M(config, ICAmodels),
        # Choose which ICAmodels M to analyse
        running_ICA.get_ICA_running(config, ICAmodels),
        plotting_ICA.get_ICA_plotting(config, ICAmodels),
        # Running the DESeq2 analysis from the config file
        DEGs.get_DESeq(config, ICAmodels),
        # Generating Supplementary Data 3 file
        refseq_noexon("results/SupplementaryData3.tsv")
