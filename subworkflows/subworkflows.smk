subworkflow refseq_noexon:
    workdir:
        "subworkflows/refseq_noexon"
    snakefile:
        "subworkflows/refseq_noexon/Snakefile"
    configfile:
        "subworkflows/refseq_noexon/config.json"


subworkflow rna_seq_cartesian_product:
    workdir:
        "subworkflows/rna_seq_cartesian_product"
    snakefile:
        "subworkflows/rna_seq_cartesian_product/Snakefile"
    configfile:
        "subworkflows/rna_seq_cartesian_product/config.json"


subworkflow rna_seq_cartesian_product_replicate:
    workdir:
        "subworkflows/rna_seq_cartesian_product_replicate"
    snakefile:
        "subworkflows/rna_seq_cartesian_product_replicate/Snakefile"
    configfile:
        "subworkflows/rna_seq_cartesian_product_replicate/config.json"
