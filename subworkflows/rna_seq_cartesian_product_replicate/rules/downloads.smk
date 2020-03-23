
def get_ebi_ftp(wildcards):
    """ Using a SRA dataset ID, generates to ftp download URL """
    data_id = wildcards.data_id
    url = ["ftp://ftp.sra.ebi.ac.uk/vol1/fastq", data_id[:6], '00' + data_id[-1], data_id, data_id]
    return '/'.join(url)


def get_fastq_R1():
    """ Generates a list of FASTQ R1 files """
    file_path = "data/datasets/{tissue}/{dataset}.R1.fastq.gz"
    dataset_dict = config['datasets']

    files = dict()
    for tissue, datasets in dataset_dict.items():
        for dataset in datasets:
            files[dataset] = file_path.format(tissue=tissue, dataset=dataset)
    return files


rule download_datasets:
    """ Download all sequencing datasets """
    output:
        R1 = config['path']['raw_fastq']['R1'],
        R2 = config['path']['raw_fastq']['R2']
    params:
        link = get_ebi_ftp
    shell:
        "wget --quiet -O {output.R1} {params.link}_1.fastq.gz && "
        "wget --quiet -O {output.R2} {params.link}_2.fastq.gz "


rule download_annotations:
    """ Download the different annotations as .GTF """
    output:
        "data/references/{annotation}.gtf"
    params:
        link = lambda wildcards:
            config["download"][wildcards.annotation]
    shell:
        "wget --quiet -O {output}.gz {params.link} && "
        "gunzip {output}.gz "


rule download_genome:
    """ Download genome """
    output:
        config['path']['genome']
    shell:
        "wget --quiet -O {output}.gz {config[download][genome]} && "
        "gunzip {output}.gz"


rule translateRefseqChr:
    """
        Transform NCBI chr ID to numeric value.
        Example : NC_000001.11 -> 1
    """
    input:
        "data/references/refseq.gtf"
    output:
        "data/references/clean_refseq_gtf.tkn"
    log:
        "logs/translateRefseqChr/gtf.log"
    script:
        "../scripts/downloads/translateRefseqChr.py"


rule extractGTFgeneID_ensembl:
    """
        From an Ensembl GTF, create a two column file with gene symbol and
        Ensembl ID
    """
    input:
        gtf = "data/references/ensembl{ens_version}.gtf"
    output:
        "data/references/ensembl{ens_version}_gene.tsv"
    shell:
        "paste "
        "<( cat {input.gtf} | "
        "awk '{{if ($3 ==\"gene\") {{print $0}}}}' FS=\"\t\" | "
        "grep -oP 'gene_name \"\K.*?(?=\";)' ) "
        "<( cat {input.gtf} | "
        "awk '{{if ($3 ==\"gene\") {{print $0}}}}' FS=\"\t\" | "
        "grep -oP 'gene_id \"\K.*?(?=\";)' ) "
        "> {output}"


rule extractGTFgeneID_refseq:
    """
        From a RefSeq GTF, create a three column file with gene symbol, RefSeq
        ID and HGNC ID.
    """
    input:
        "data/references/clean_refseq_gtf.tkn",
        gtf = "data/references/refseq.gtf"
    output:
        refseq_gene = "data/references/refseq_gene.tsv"
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/downloads/extractGTFgeneID_refseq.py"


rule download_HGNC:
    """
        Downloads an HGNC table with the custom download tool. Downloads data
        provided by HGNC and by external sources.
    """
    output:
        temp("data/references/hgnc.txt.temp")
    shell:
        "wget --quiet -O {output} \"{config[download][HGNC]}\""


rule clean_HGNC:
    """
        For some HGNC ID, data from HGNC and external sources is not equal. If
        one entry for a source is missing, the other source is used. If the two
        sources do not have the same data, HGNC is prioritized.
    """
    input:
        hgnc = "data/references/hgnc.txt.temp"
    output:
        hgnc = "data/references/hgnc.txt"
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/downloads/clean_HGNC.py"


rule get_datasets_depth:
    """
        Calculates a scaling metric based on the number of lines in the FASTQ
        files to approximate depth of the dataset.
    """
    input:
        **get_fastq_R1()
    output:
        datasets_depth = "data/datasets/datasets_depth.tsv"
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/downloads/get_datasets_depth.py"


rule indexGenome:
    """ Index the genome using Samtools """
    input:
        config['path']['genome']
    output:
        config['path']['genome'] + '.fai'
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools faidx {input}"
