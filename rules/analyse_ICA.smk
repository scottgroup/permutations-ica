
rule GO_analysis:
    input:
        gene_list = "results/{ICA_path}/gene_list/comp_sigma{sigma}/.tkn"
    output:
        plot = "results/{ICA_path}/gene_list/comp_sigma{sigma}/GO_comp{comp}.svg",
    params:
        path = "results/{ICA_path}/gene_list/comp_sigma{sigma}/comp_{comp}.txt"
    conda:
        "../envs/gprofiler.yaml"
    script:
        "../scripts/analyse_ICA/GO_analysis.py"
