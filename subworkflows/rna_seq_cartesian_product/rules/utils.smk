import get_inputs


rule combine_quantification:
    input:
        get_inputs.quantification_results(config, '{annotation}'),
        gene_list = "data/references/{annotation}_gene.tsv"
    output:
        out_file = "results/cartesian_product/raw_noNaN_{annotation}.tsv",
        out_NaN_file = "results/cartesian_product/raw_NaN_{annotation}.tsv"
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/utils/combine_quantification.py"


rule merge_quantifications:
    input:
        quants = expand(
            "results/cartesian_product/raw_{{isNaN}}_{annotation}.tsv",
            annotation=config['tools']['annotation']
        ),
        hgnc = "data/references/hgnc.txt"
    output:
        out_file = "results/cartesian_product/tissues_{isNaN}.tsv"
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/utils/merge_quantifications.py"
