

def get_blat_results(wildcards):
    """ """
    split_ids, = glob_wildcards("data/fasta/processed_pseudogene.{id}.fa")
    if len(split_ids) > 0:
        return expand(rules.run_blat.output.psl, id=split_ids)
    else:
        return []


rule make_ooc_file:
    input:
        transcriptome = "data/fasta/protein_coding.fa"
    output:
        ooc = "data/protein_coding.11.ooc"
    conda:
        "../envs/blat.yaml"
    shell:
        "blat {input.transcriptome} /dev/null /dev/null -makeOoc={output.ooc}"


rule run_blat:
    input:
        rules.splitting_queryfasta.output.tkn,
        ooc = rules.make_ooc_file.output.ooc,
        transcriptome = "data/fasta/protein_coding.fa",
    output:
        psl = "results/blat/processed_pseudogene.{id}.blast8"
    conda:
        "../envs/blat.yaml"
    params:
        query = "data/fasta/processed_pseudogene.{id}.fa",
    shell:
        "blat {input.transcriptome} {params.query} "
        "-fine -noHead -out=blast8 "
        "{output.psl}"


rule compile_blat:
    input:
        blats = get_blat_results,
        gene_map = rules.extract_gene_id.output.gene_map,
        tkn = rules.splitting_queryfasta.output.tkn
    output:
        score = "results/blat_score.tsv"
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/compile_blat.py"
