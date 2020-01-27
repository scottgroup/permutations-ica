subworkflow refseq_noexon:
    workdir:
        "subworkflows/refseq_noexon"
    snakefile:
        "subworkflows/refseq_noexon/Snakefile"
    configfile:
        "subworkflows/refseq_noexon/config.json"


subworkflow pseudogene_parent:
    workdir:
        "subworkflows/pseudogene_parent"
    snakefile:
        "subworkflows/pseudogene_parent/Snakefile"
    configfile:
        "subworkflows/pseudogene_parent/config.json"


subworkflow rna_seq_cartesian_product:
    workdir:
        "subworkflows/rna_seq_cartesian_product"
    snakefile:
        "subworkflows/rna_seq_cartesian_product/Snakefile"
    configfile:
        "subworkflows/rna_seq_cartesian_product/config.json"
