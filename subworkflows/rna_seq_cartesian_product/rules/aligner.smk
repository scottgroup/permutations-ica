
rule bowtie2_index:
    input:
        genome = config['path']['genome']
    output:
        bowtie2_index = directory("data/indexes/bowtie2"),
    params:
        bowtie2_index = "data/indexes/bowtie2/bowtie2"
    conda:
        "../envs/tophat2.yaml"
    threads:
        32
    log:
        "logs/tophat2/bowtie2_index.log"
    shell:
        "mkdir -p {output.bowtie2_index} && "
        "bowtie2-build "
        "--threads {threads} "
        "{input.genome} "
        "{params.bowtie2_index} &> {log}"


rule tophat2_align:
    input:
        fastq1 = "results/rnaseq/{trimmer}/{datasets}/{rep}/R1.fastq",
        fastq2 = "results/rnaseq/{trimmer}/{datasets}/{rep}/R2.fastq",
        bowtie2_index = "data/indexes/bowtie2"
    output:
        "results/rnaseq/tophat2/{datasets}/{rep}/{trimmer}/aligned.bam"
    params:
        bowtie2_index = "data/indexes/bowtie2/bowtie2",
        out_dir = "results/rnaseq/tophat2/{datasets}/{rep}/{trimmer}"
    conda:
        "../envs/tophat2.yaml"
    threads:
        32
    log:
        "logs/tophat2/align_{datasets}_{rep}_{trimmer}.log"
    shell:
        "tophat "
        "--num-threads {threads} "
        "--output-dir {params.out_dir} "
        "--library-type fr-unstranded "
        "{params.bowtie2_index} "
        "{input.fastq1} {input.fastq2} &> {log} && "
        "mv {params.out_dir}/accepted_hits.bam {output} && "
        "rm {params.out_dir}/unmapped.bam"


rule HISAT2_index:
    input:
        genome = config['path']['genome']
    output:
        directory("data/indexes/HISAT2/")
    conda:
        "../envs/HISAT2.yaml"
    threads:
        32
    log:
        "logs/HISAT2/HISAT2_index.log"
    shell:
        "hisat2-build "
        "--threads {threads} "
        "{input.genome} "
        "{output} &> {log}"


rule HISAT2_align:
    input:
        fastq1 = "results/rnaseq/{trimmer}/{datasets}/{rep}/R1.fastq",
        fastq2 = "results/rnaseq/{trimmer}/{datasets}/{rep}/R2.fastq",
        idx = "data/indexes/HISAT2/"
    output:
        "results/rnaseq/HISAT2/{datasets}/{rep}/{trimmer}/aligned.bam"
    params:
        out_dir = "results/rnaseq/HISAT2/{datasets}/{rep}/{trimmer}"
    conda:
        "../envs/HISAT2.yaml"
    threads:
        32
    log:
        "logs/HISAT2/align_{datasets}_{rep}_{trimmer}.log"
    shell:
        "mkdir -p {params.out_dir} && "
        "(hisat2 "
        "--threads {threads} "
        "-x {input.idx} -1 {input.fastq1} -2 {input.fastq2} | "
        "samtools view -Sbh -o {output} ) &> {log} "


rule STAR_index:
    input:
        "data/references/clean_refseq_gtf.tkn",
        genome = config['path']['genome'],
        gtf = "data/references/{annotation}.gtf"
    output:
        directory("data/indexes/STAR/{annotation}/")
    conda:
        "../envs/STAR.yaml"
    threads:
        32
    params:
        sjdbOverhang = "99"
    log:
        "logs/STAR/STAR_generateGenome_{annotation}.log"
    shell:
        "STAR "
        "--runMode genomeGenerate "
        "--runThreadN {threads} "
        "--genomeDir {output} "
        "--outFileNamePrefix {output} "
        "--genomeFastaFiles {input.genome} "
        "--sjdbGTFfile {input.gtf} "
        "--sjdbOverhang {params.sjdbOverhang} "
        "&> {log}"


rule STAR_align:
    input:
        fastq1 = "results/rnaseq/{trimmer}/{datasets}/{rep}/R1.fastq",
        fastq2 = "results/rnaseq/{trimmer}/{datasets}/{rep}/R2.fastq",
        genomeDir = "data/indexes/STAR/{annotation}/"
    output:
        "results/rnaseq/STAR/{datasets}/{rep}/{trimmer}/{annotation}/aligned.bam"
    params:
        out_dir = "results/rnaseq/STAR/{datasets}/{rep}/{trimmer}/{annotation}",
        sjdbOverhang = "99"
    conda:
        "../envs/STAR.yaml"
    threads:
        32
    log:
        "logs/STAR/align_{datasets}_{rep}_{trimmer}_{annotation}.log"
    shell:
        "mkdir -p {params.out_dir} && "
        "STAR "
        "--runThreadN {threads} "
        "--genomeDir {input.genomeDir} "
        "--readFilesIn {input.fastq1} {input.fastq2} "
        "--sjdbOverhang {params.sjdbOverhang} "
        "--outFileNamePrefix {params.out_dir}/ "
        "--outSAMtype BAM Unsorted "
        "--outStd Log &> {log} && "
        "mv {params.out_dir}/Aligned.out.bam {output}"
