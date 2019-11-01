
rule download_annotation:
    output:
        gtf = "data/references/ensembl.gtf.gz"
    params:
        link = config['path']['gtf']
    shell:
        "wget --quiet -O {output.gtf} {params.link}"


rule download_genome:
    output:
        genome = "data/references/genome.fa"
    params:
        link = config['path']['genome']
    shell:
        "wget --quiet -O {output.genome}.gz {params.link} && "
        "gunzip {output}.gz "
