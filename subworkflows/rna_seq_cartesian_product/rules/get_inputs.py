from snakemake.io import expand

def download_list(config):
    """
    Listing all downloads.
    Useful to run independently if cluster nodes have not internet access.
    """
    file_list = list()
    # Genome
    file_list.append("data/references/genome.fa")
    # HGNC file
    file_list.append("data/references/hgnc.txt")
    # Annotations
    file_list.extend(
        expand(
            "data/references/{annotation}.gtf",
            annotation=config['tools']['annotation']
        )
    )
    # Correcting refseq
    file_list.append("data/references/clean_refseq_gtf}.tkn")

    # For each tissue
    for tissue in config['datasets'].keys():
        file_list.extend(
            expand(
                "data/datasets/{tissue}/{data_id}.R1.fastq.gz",
                tissue=tissue,
                data_id=config['datasets'][tissue]
            )
        )
    return file_list


def quantification_results(config, use_annotation=False):
    """
    List all quantification results. Marks the end of the RNA-seq pipeline
    """
    if not use_annotation:
        use_annotation = config['tools']['annotation']

    file_list = list()
    # For each tissue
    for tissue in config['datasets'].keys():
        _counts = expand(
            "results/rnaseq/{quantifier}/{datasets}/{data_id}/{trimmer}/{aligner}/{annotation}/counts.tsv",
            datasets=tissue,
            data_id=config['datasets'][tissue],
            trimmer=config['tools']['trimmer'],
            annotation=use_annotation,
            aligner=config['tools']['aligner'],
            quantifier=config['tools']['quantifier']
        )
        file_list.extend(_counts)
    return file_list
