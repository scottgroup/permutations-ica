
def get_sorted_bam(wildcards):
    """ Format BAM file respective to the aligner """
    if wildcards["aligner"] == 'STAR':
        st = "results/rnaseq/{aligner}/{datasets}/{rep}/{trimmer}/{annotation}/aligned.sorted.bam"
    else:
        st = "results/rnaseq/{aligner}/{datasets}/{rep}/{trimmer}/aligned.sorted.bam"
    return st.format(**wildcards)


def get_all_sorted_bam(wildcards):
    """ Generate all BAM files for the different software in the study """
    sorted_bam = list()
    for tissue, datasets in config['datasets'].items():
        wildcard = {"datasets": tissue}
        for dataset in datasets:
            wildcard["rep"] = dataset
            for aligner in config['tools']['aligner']:
                wildcard["aligner"] = aligner
                for trimmer in config['tools']['trimmer']:
                    wildcard["trimmer"] = trimmer
                    for annotation in config['tools']['annotation']:
                        wildcard["annotation"] = annotation
                        sorted_bam.append(get_sorted_bam(wildcard))
    return list(set(sorted_bam))


def get_all_sorted_bai(wildcards):
    """ For all BAM files, also create a .BAM.BAI file """
    return [file + '.bai' for file in get_all_sorted_bam(wildcards)]


rule sortBAM:
    """ Sorts the BAM files """
    input:
        "results/rnaseq/{aligner}/{datasets}/{rep}/{trimmer_annot}/aligned.bam"
    output:
        "results/rnaseq/{aligner}/{datasets}/{rep}/{trimmer_annot}/aligned.sorted.bam"
    conda:
        "../envs/samtools.yaml"
    threads:
        32
    log:
        "logs/samtools_sort/{datasets}_{rep}_{trimmer_annot}_{aligner}.log"
    shell:
        "samtools sort -@ {threads} -o {output} {input} && "
        "rm {input}"


rule indexBAM:
    """ Indexes the sorted BAM files """
    input:
        "results/rnaseq/{aligner}/{datasets}/{rep}/{trimmer_annot}/aligned.sorted.bam"
    output:
        "results/rnaseq/{aligner}/{datasets}/{rep}/{trimmer_annot}/aligned.sorted.bam.bai"
    conda:
        "../envs/samtools.yaml"
    threads:
        1
    log:
        "logs/samtools_index/{datasets}_{rep}_{trimmer_annot}_{aligner}.log"
    shell:
        "samtools index {input}"


rule cufflinks:
    """ Quantifies the alignment through Cufflinks """
    input:
        gtf = "data/references/{annotation}.gtf",
        bam = get_sorted_bam
    output:
        "results/rnaseq/cufflinks/{datasets}/{rep}/{trimmer}/{aligner}/{annotation}/counts.tsv"
    params:
        out_dir = "results/rnaseq/cufflinks/{datasets}/{rep}/{trimmer}/{aligner}/{annotation}"
    conda:
        "../envs/cufflinks.yaml"
    threads:
        32
    log:
        "logs/cufflinks/{datasets}_{rep}_{trimmer}_{aligner}_{annotation}.log"
    shell:
        "cuffquant "
        "--num-threads {threads} "
        "--output-dir {params.out_dir} "
        "--library-type fr-unstranded "
        "{input.gtf} "
        "{input.bam} &> {log} "
        "&& cuffnorm "
        "--num-threads {threads} "
        "--output-dir {params.out_dir} "
        "--library-type fr-unstranded "
        "{input.gtf} "
        "{params.out_dir}/abundances.cxb {params.out_dir}/abundances.cxb "
        " &>> {log} && "
        "sed '1d' {params.out_dir}/genes.count_table | "
        "awk '{{print $1,$2}}' > {output} "
        " && rm {params.out_dir}/abundances.cxb "
        "{params.out_dir}/cds.* {params.out_dir}/genes.* "
        "{params.out_dir}/isoforms.* {params.out_dir}/tss_groups.* "
        "{params.out_dir}/run.info {params.out_dir}/samples.table"


rule htseq:
    """ Quantifies the alignment through HTSeq """
    input:
        gtf = "data/references/{annotation}.gtf",
        bam = get_sorted_bam
    output:
        "results/rnaseq/htseq/{datasets}/{rep}/{trimmer}/{aligner}/{annotation}/counts.tsv"
    params:
        out_dir = "results/rnaseq/htseq/{datasets}/{rep}/{trimmer}/{aligner}/{annotation}"
    conda:
        "../envs/htseq.yaml"
    threads:
        32
    log:
        "logs/htseq/{datasets}_{rep}_{trimmer}_{aligner}_{annotation}.log"
    shell:
        "htseq-count "
        "--format bam "
        "--mode intersection-nonempty "
        "--order pos "
        "--type exon "
        "--idattr gene_id "
        "--stranded=no "
        "{input.bam} "
        "{input.gtf} > {output} 2> {log}"


rule featureCounts:
    """ Quantifies the alignment through featureCounts """
    input:
        gtf = "data/references/{annotation}.gtf",
        bam = get_sorted_bam
    output:
        "results/rnaseq/featureCounts/{datasets}/{rep}/{trimmer}/{aligner}/{annotation}/counts.tsv"
    conda:
        "../envs/subread.yaml"
    threads:
        32
    log:
        "logs/featureCounts/{datasets}_{rep}_{trimmer}_{aligner}_{annotation}.log"
    shell:
        "featureCounts "
        "-T {threads} "
        "-a {input.gtf} "
        "-o {output}.temp "
        "{input.bam} &> {log} && "
        "sed '1d' {output}.temp | awk '{{print $1,$6,$7}}' > {output} && "
        "mv {output}.temp.summary {output}.summary && "
        "rm {output}.temp"


rule getGeneCoverage:
    """
        For a specific gene, using the ENGS00000XXXXX gene ID, extracts read
        coverage from the different pipelines, using Ensembl 98 gene
        coordinates.
    """
    input:
        bams = get_all_sorted_bam,
        bais = get_all_sorted_bai,
        datasets_depth = rules.get_datasets_depth.output.datasets_depth,
        gtf = "data/references/ensembl98.gtf"
    output:
        gene_coverage = "results/rnaseq/geneCoverage/{gene}.tsv",
        gene_BAM_tkn = "results/geneCoverage/{gene}.tkn"
    conda:
        "../envs/python.yaml"
    threads:
        1
    script:
        "../scripts/quantification/getGeneCoverage.py"
