
rule GO_analysis:
    input:
        gene_list = "results/{ICA_path}/gene_list/comp_sigma{sigma}/.tkn"
    output:
        plot_up = "results/{ICA_path}/gene_list/comp_sigma{sigma}/GO_comp{comp}_up.png",
        plot_down = "results/{ICA_path}/gene_list/comp_sigma{sigma}/GO_comp{comp}_down.png"
    params:
        path = "results/{ICA_path}/gene_list/comp_sigma{sigma}/comp_{comp}.txt"
    conda:
        "../envs/gprofiler.yaml"
    script:
        "../scripts/analyse_ICA/GO_analysis.py"
