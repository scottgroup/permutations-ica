

rule generating_pseudogene_list:
    input:
        gtf = rules.download_annotation.output.gtf
    output:
        pseudo_list = "data/{biotype}.txt"
    params:
        bio = lambda w: "\|".join(config['biotype']['{biotype}'.format(**w)])
    shell:
        "zcat {input.gtf} | "
        "awk '$3==\"gene\" {{print $0}}' | "
        "grep -w '{params.bio}' | "
        "grep -oP 'gene_id \"\K.*?(?=\")' "
        "> {output.pseudo_list}"


rule generating_sub_gtf:
    input:
        gtf = rules.download_annotation.output.gtf,
        list = "data/{biotype}.txt"
    output:
        sub_gtf = "data/{biotype}.gtf"
    shell:
        "zcat {input.gtf} | "
        "grep -f {input.list} "
        "> {output.sub_gtf}"


rule generating_fasta:
    input:
        genome = rules.download_genome.output.genome,
        gtf = rules.generating_sub_gtf.output.sub_gtf
    output:
        transcriptome = "data/fasta/{biotype}.fa"
    conda:
        "../envs/gffread.yaml"
    shell:
        "gffread -w {output.transcriptome} -g {input.genome} {input.gtf}"


checkpoint splitting_queryfasta:
    input:
        fasta = "data/fasta/processed_pseudogene.fa"
    output:
        tkn = "data/fasta/processed_pseudogene.tkn"
    conda:
        "../envs/pyfasta.yaml"
    params:
        file = "%(fasta)s/%(seqid)s.fa"
    shell:
        "pyfasta split -n 200 {input.fasta} && "
        "touch {output.tkn}"


rule extract_gene_id:
    input:
        gtf = rules.download_annotation.output.gtf
    output:
        gene_map = "data/references/gene_id.tsv"
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/extract_gene_id.py"
