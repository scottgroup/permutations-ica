import get_inputs


rule combine_quantification:
    """
        Generates a files with all quantification from one annotation for the
        different pipelines.

        In the NaN file, to be reported, genes must have a quantification for
        each pipeline. In the noNaN, if a gene has quantification in at least
        one pipeline, it will be reported, and will have a 0 value for each
        other pipeline where it is not present.
    """
    input:
        get_inputs.quantification_results(config, '{strand}', '{annotation}'),
        gene_list = "data/references/{annotation}_gene.tsv"
    output:
        out_file = "results/cartesian_product/raw_{strand}_noNaN_{annotation}.tsv",
        out_NaN_file = "results/cartesian_product/raw_{strand}_NaN_{annotation}.tsv"
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/utils/combine_quantification.py"


rule merge_quantifications:
    """
        Merge the different quantifications from rules.combine_quantification
        using HGNC to bridge between IDs.
    """
    input:
        quants = expand(
            "results/cartesian_product/raw_{{strand}}_{{isNaN}}_{annotation}.tsv",
            annotation=config['tools']['annotation']
        ),
        hgnc = "data/references/hgnc.txt"
    output:
        out_file = "results/cartesian_product/tissues_{strand}_{isNaN}.tsv"
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/utils/merge_quantifications.py"
