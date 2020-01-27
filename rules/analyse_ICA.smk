from functools import partial

wildcard_constraints:
    side = "(up|down|both)",


def getQuantFolder(wildcards):
    if wildcards["aligner"] == 'STAR':
        st = "subworkflows/rna_seq_cartesian_product/results/rnaseq/{aligner}/{datasets}/{rep}/{trimmer}/{annotation}/"
    else:
        st = "subworkflows/rna_seq_cartesian_product/results/rnaseq/{aligner}/{datasets}/{rep}/{trimmer}/"
    return st.format(**wildcards)


def getAllQuantFolder(wildcards):
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
                        sorted_bam.append(getQuantFolder(wildcard))
    return list(set(sorted_bam))


def getGeneFiles(path, wildcards):
    """ """
    gene_list = "results/{ICA_path}/sigma_{sigma}/gene_list/comp_{comp}.txt"
    genes = list()
    with open(gene_list.format(**wildcards), 'r') as f:
        for line in f.readlines():
            if line.strip() == ">Positive genes":
                up = True
            elif line.strip() == ">Negative genes":
                up = False
            elif wildcards.side == 'up' and up:
                genes.append(line.strip())
            elif wildcards.side == 'down' and not up:
                genes.append(line.strip())

    files = list()
    for gene in genes:
        files.append(rna_seq_cartesian_product(path.format(gene=gene)))
    return files


rule GO_analysis:
    """
        For a component, generates GO analysis and gene weights distribution.

        ori -> [up, down, both] specifies which analysis to include. Up and down
            are always reported independently, on their respective side of the
            gene weights distribution plot.
    """
    input:
        gene_list = "results/{ICA_path}/sigma_{sigma}/gene_list/.tkn",
        components_mean = "results/{ICA_path}/components_mean.tsv"
    output:
        plot = "results/{ICA_path}/sigma_{sigma}/gene_list/GO_comp{comp}_{ori}.svg",
    params:
        path = "results/{ICA_path}/sigma_{sigma}/gene_list/comp_{comp}.txt"
    conda:
        "../envs/gprofiler.yaml"
    script:
        "../scripts/analyse_ICA/GO_analysis.py"


rule getMetaGene:
    """
        For every gene in a component side, retrieve the previously calculated
        gene coverage and compute the average profile of the various acceptor
        and donor sites. Parameters for acceptable sites and length of coverage
        are described in the script. The profiles are separated by aligners.

        side -> [up, down]
    """
    input:
        genes = partial(getGeneFiles, "results/rnaseq/geneCoverage/{gene}.tsv"),
        gtf = rna_seq_cartesian_product("data/references/ensembl98.gtf")
    output:
        mean = "results/{ICA_path}/metaGene/comp_{comp}_{side}_mean.tsv",
        std = "results/{ICA_path}/metaGene/comp_{comp}_{side}_std.tsv"
    conda:
        "../envs/ICA_python.yaml"
    script:
        "../scripts/analyse_ICA/getMetaGene.py"


rule plotMetaGene:
    """
        Plot a metaGene for specified compononent and side.
    """
    input:
        mean = rules.getMetaGene.output.mean,
        std = rules.getMetaGene.output.std
    output:
        plot = "results/{ICA_path}/metaGene/comp_{comp}_{side}.svg"
    conda:
        "../envs/ICA_python.yaml"
    script:
        "../scripts/analyse_ICA/plotMetaGene.py"


rule getFlagstatOtherChr:
    """
        For all gene in a component side, calculate the % of mate mapped to
        other chr using Flagstat. Quantifies and store results to be plotted by
        plotFlagstatOtherChr rule.
    """
    input:
        gene_BAM_tkn = partial(getGeneFiles, "results/geneCoverage/{gene}.tkn"),
        folders = getAllQuantFolder
    output:
        file = "results/{ICA_path}/metaGene/comp_{comp}_{side}_sigma{sigma}_chr.tsv"
    conda:
        "../envs/samtools.yaml"
    script:
        "../scripts/analyse_ICA/getFlagstatOtherChr.py"


rule plotFlagstatOtherChr:
    """
        Plotting boxplot of the % of mate mapped to other chr.
    """
    input:
        file = rules.getFlagstatOtherChr.output.file
    output:
        plot = "results/{ICA_path}/metaGene/comp_{comp}_{side}_sigma{sigma}_chr.svg"
    conda:
        "../envs/ICA_python.yaml"
    script:
        "../scripts/analyse_ICA/plotFlagstatOtherChr.py"


rule getGeneSeq:
    """
        Extract gene sequence using a gene feature chr, start and end columns.
        Strand is not taken into consideration, sequence from the sense strand
        will always be reported.
    """
    input:
        genome = rna_seq_cartesian_product("data/references/genome.fa"),
        genome_index = rna_seq_cartesian_product("data/references/genome.fa.fai"),
        gtf = rna_seq_cartesian_product("data/references/ensembl98.gtf")
    output:
        gene = "results/ICA/metaGenes/seqs/{gene}.fa"
    conda:
        "../envs/samtools.yaml"
    script:
        "../scripts/analyse_ICA/getGeneSeq.py"


rule plotGeneVsPseudogene:
    """
        Generates a plot displaying a gene, its pseudogene and their alignment
        profiles from the diverse pipelines included in the study.

        A transcript must be specify for the main gene.
        Pseudogenes are expected to have only one transcript.
    """
    input:
        gene = "subworkflows/rna_seq_cartesian_product/results/rnaseq/geneCoverage/{gene}.tsv",
        pseudo = "subworkflows/rna_seq_cartesian_product/results/rnaseq/geneCoverage/{pseudogene}.tsv",
        gene_seq = "results/ICA/metaGenes/seqs/{gene}.fa",
        pseudo_seq = "results/ICA/metaGenes/seqs/{pseudogene}.fa",
        gtf = rna_seq_cartesian_product("data/references/ensembl98.gtf")
    output:
        plot = "results/ICA/metaGenes/{gene}_vs_{pseudogene}_{transcript_id}.svg"
    conda:
        "../envs/ICA_python.yaml"
    script:
        "../scripts/analyse_ICA/plotGeneVsPseudogene.py"
