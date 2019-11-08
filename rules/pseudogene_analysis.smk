
rule download_pseudogene_parents:
    output:
        file = "data/references/pseudogene_parents.txt"
    params:
        link = lambda w: "{pseudogenes_parents}".format(**config['url'])
    shell:
        "wget --quiet -O {output.file} {params.link}"


rule download_pseudogene_biotype:
    output:
        file = "data/references/pseudogene_biotype.txt"
    params:
        link = lambda w: "{pseudogenes_biotype}".format(**config['url'])
    shell:
        "wget --quiet -O {output.file} {params.link}"


rule inspect_pseudogenes_parents:
    input:
        parent = rules.download_pseudogene_parents.output.file,
        biotype = rules.download_pseudogene_biotype.output.file,
        gene_list_tkn = "results/ICA/{ICAmethod}/{dataset}/{ICA_run}/gene_list/comp_sigma{sigma}/.tkn",
        data = lambda w: config['ICA_datasets']['{dataset}'.format(**w)]['params']['counts'],
        stuff=pseudogene_parent("results/blat_score.tsv")
    output:
        "results/ICA/{ICAmethod}/{dataset}/{ICA_run}/sigma_{sigma}/pseudogene_comp{comp}.png"
    params:
        gene_list = "results/ICA/{ICAmethod}/{dataset}/{ICA_run}/gene_list/comp_sigma{sigma}/comp_{comp}.txt"
    conda:
        "../envs/ICA_python.yaml"
    script:
        "../scripts/pseudogene_analysis/inspect.py"
