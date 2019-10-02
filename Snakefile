configfile: "config.json"

import rules.running_ICA as running_ICA
import rules.plotting_ICA as plotting_ICA

# Defining wildcards constraints
wildcard_constraints:
    dataset = "({})".format("|".join(config["ICA_datasets"].keys())),
    ICAmethod = "sklearnFastICA",
    std="[0-9]+"


# Including rules
include: "rules/running_ICA.smk"
include: "rules/plotting_ICA.smk"

# Defining model to run
ICAruns = [
    "counts_NaN", "counts_noNaN", 
    "counts_NaN_dualAnnot", "counts_noNaN_dualAnnot",
    "counts_NaN_dualAnnot_heart", "counts_noNaN_dualAnnot_heart"
]


rule all:
    input:
        running_ICA.get_ICA_running(config, ICAruns),
        plotting_ICA.get_ICA_plotting(config, ICAruns)
