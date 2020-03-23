
def get_cutadapt_params():
    """ Formatting trimming parameters for Cutadapt """
    return "-q {phred} --minimum-length {minlen}".format(
        phred=config["trimmer"]["min_phred"],
        minlen=config["trimmer"]["min_length"]
    )


def get_trimmomatic_params():
    """ Formatting trimming parameters for Trimmomatic """
    return [
        "SLIDINGWINDOW:5:{phred}".format(phred=config["trimmer"]["min_phred"]),
        "MINLEN:{min_len}".format(min_len=config["trimmer"]["min_length"])
    ]


rule cutadapt:
    """ Trimming using Cutadapt """
    input:
        fastq1 = config['path']['raw_fastq']['R1'],
        fastq2 = config['path']['raw_fastq']['R2']
    output:
        fastq1 = "results/rnaseq/cutadapt/{datasets}/{data_id}/R1.fastq",
        fastq2 = "results/rnaseq/cutadapt/{datasets}/{data_id}/R2.fastq"
    threads:
        32
    params:
        get_cutadapt_params()
    conda:
        "../envs/cutadapt.yaml"
    log:
        "logs/cutadapt/{datasets}_{data_id}.log"
    shell:
        "cutadapt "
        "--cores={threads} "
        "{params} "
        "-o {output.fastq1} "
        "-p {output.fastq2} "
        "{input.fastq1} {input.fastq2} "
        " &> {log}"


rule trimmomatic:
    """ Trimming using Trimmomatic """
    input:
        r1 = config['path']['raw_fastq']['R1'],
        r2 = config['path']['raw_fastq']['R2']
    output:
        r1 = "results/rnaseq/trimmomatic/{datasets}/{data_id}/R1.fastq",
        r2 = "results/rnaseq/trimmomatic/{datasets}/{data_id}/R2.fastq",
        r1_unpaired = "results/rnaseq/trimmomatic/{datasets}/{data_id}/unpairedR1.fastq",
        r2_unpaired = "results/rnaseq/trimmomatic/{datasets}/{data_id}/unpairedR2.fastq",
    params:
        trimmer = get_trimmomatic_params()
    conda:
        "../envs/trimmomatic.yaml"
    threads:
        32
    log:
        "logs/trimmomatic/{datasets}_{data_id}.log"
    shell:
        "trimmomatic PE "
        "-threads {threads} "
        "{input.r1} {input.r2} "
        "{output.r1} {output.r1_unpaired} "
        "{output.r2} {output.r2_unpaired} "
        "{params.trimmer} "
        "&> {log}"
