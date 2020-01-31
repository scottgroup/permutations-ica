from functools import partial

from genomic_overlap_analysis import getGeneOverlaps
from genomic_overlap_analysis import getAdjacentGenes

wildcard_constraints:
    annotation = "({})".format("|".join(config["tools"]["annotation"])),
    feature = "(gene|exon)",
    comp = "[0-9]+",
    HGNCgene="HGNC:[0-9]+"


rule get_local_gtf:
    """ Move and zip a copy of the GTF file for each annotation """
    input:
        gtf = rna_seq_cartesian_product("data/references/{annotation}.gtf")
    output:
        gtf = "data/references/{annotation}.gtf.gz"
    shell:
        "cat {input.gtf} | gzip -9 > {output.gtf}"


rule translateRefseqChr:
    """ Change RefSeq chr names to match Ensembl chr names """
    input:
        "data/references/refseq.gtf.gz"
    output:
        "data/references/.tkn_cleanrefseq"
    script:
        "../scripts/genomic_overlap/translateRefseq.py"


rule getOnlyFeature:
    """ Extract rows from GTF that match the feature type """
    input:
        gtf = rules.get_local_gtf.output.gtf,
        refseq_tkn = "data/references/.tkn_cleanrefseq"
    output:
        "results/genomic_overlap/annotation_{feature}/{annotation}.gtf.gz"
    shell:
        "set +o pipefail; "
        "zcat {input.gtf} | head -n 50 | grep '#!' | gzip -9 > {output}; "
        "zcat {input.gtf} | "
        "awk 'BEGIN {{FS=\"\t\"}}; {{ if ($3 == \"{wildcards.feature}\") {{print $0}}}}' | "
        "gzip -9 > {output}"


rule getBed6_ensembl:
    """ Generate BED6 from GTF for Ensembl """
    input:
        gene_gtf = "results/genomic_overlap/annotation_{feature}/ensembl{ens_version}.gtf.gz"
    output:
        bed = "results/genomic_overlap/bed6/ensembl{ens_version}_{feature}.bed"
    shell:
        "zcat {input.gene_gtf} | "
        "awk '{{print $1,$4-1,$5,$9,1,$7}}' OFS=\"\t\" FS=\"[\t;]\" | "
        "tr -d \"gene_id \" | tr -d '\"' > {output.bed}"


rule getBed6_refseq:
    """ Generate BED6 from GTF for RefSeq """
    input:
        gene_gtf = "results/genomic_overlap/annotation_{feature}/refseq.gtf.gz"
    output:
        bed = "results/genomic_overlap/bed6/refseq_{feature}.bed"
    shell:
        "zcat {input.gene_gtf} | "
        "awk 'BEGIN {{FS=\"\t\"}} match($9, /GeneID:[0-9]*/) "
        "{{print $1,$4-1,$5,\"GeneID:\"substr($9, RSTART+7, RLENGTH-7),1,$7}}'"
        " OFS=\"\t\" > {output.bed}"


rule getBedtoolsIntersect:
    """
        Calculate the bedtools intersection of a specific feature in two
        different annotations.
    """
    input:
        bed1 = "results/genomic_overlap/bed6/{annot1}_{feature}.bed",
        bed2 = "results/genomic_overlap/bed6/{annot2}_{feature}.bed"
    output:
        temp("results/genomic_overlap/{feature}_{annot1}_{annot2}_intersect.tsv")
    conda:
        "../envs/bedtools.yaml"
    shell:
        "bedtools intersect -a {input.bed1} -b {input.bed2} -wao | "
        "awk '{{if ($7 != \".\") {{print $0}}}}' > {output}"


rule extractOverlap:
    """
        Transforms the intersect files in overlaps files.

        column 1 -> gene1
        column 2 -> gene2
        column 3 -> % of gene1 overlapped by gene2
        column 4 -> % of gene2 overlapped by gene1
    """
    input:
        intersect = "results/genomic_overlap/{intersect}_intersect.tsv"
    output:
        overlap = "results/genomic_overlap/{intersect}_overlaps.tsv"
    shell:
        "cat {input.intersect} | "
        "awk '{{ print $4,$10,$13/($3-$2),$13/($9-$8) }}' FS='\t' OFS='\t' "
        "> {output.overlap}"


rule getGeneOverlap:
    """
        For a gene in an annotation, generate BED and overlap information
    """
    input:
        gtf = rules.get_local_gtf.output.gtf,
        refseq_tkn = "data/references/.tkn_cleanrefseq",
        exon_bed = "results/genomic_overlap/bed6/{annotation}_exon.bed"
    output:
        overlap = "results/genomic_overlap/gene_overlap/{gene}_{annotation}.tsv",
        bed = "results/genomic_overlap/gene_bed/{gene}_{annotation}.bed",
        merged_bed = "results/genomic_overlap/gene_bed_merged/{gene}_{annotation}.bed"
    params:
        dir = "results/genomic_overlap/gene/"
    conda:
        "../envs/bedtools.yaml"
    script:
        "../scripts/genomic_overlap/getGeneOverlap.py"


rule quantifyingGeneOverlap:
    """
        For a component, generates three violinplots for each significant gene
        groups.

        Violonplots show, separated by genome annotation, the percentage of exon
        overlap by other genes.
    """
    input:
        gene_list = "results/{ICA_path}/sigma_{sigma}/gene_list/comp_{comp}.txt",
        gene_list_tkn = "results/{ICA_path}/sigma_{sigma}/gene_list/.tkn",
        HGNC = rna_seq_cartesian_product("data/references/hgnc.txt"),
        overlaps = getGeneOverlaps
    output:
        plot = "results/{ICA_path}/sigma_{sigma}/comp_{comp}_overlap.svg"
    conda:
        "../envs/ICA_python.yaml"
    script:
        "../scripts/genomic_overlap/quantifyingGeneOverlap.py"


rule quantifyingUniqueGeneLength:
    """
        For a component, generates a stacked bar chart for each significant gene
        groups, separated by pairwise comparison of genome annotation.

        For each gene in each group, we extracted the gene coordinates and
        quantified the percentage of unique sequence for each annotation in a
        pairwise manner. Colored regions represent the average proportion of
        unique sequence.
    """
    input:
        gene_list = "results/{ICA_path}/sigma_{sigma}/gene_list/.tkn",
        HGNC = rna_seq_cartesian_product("data/references/hgnc.txt"),
        annotations = [
            "results/genomic_overlap/gene_ensembl92_ensembl98_overlaps.tsv",
            "results/genomic_overlap/gene_ensembl92_refseq_overlaps.tsv",
            "results/genomic_overlap/gene_ensembl98_refseq_overlaps.tsv",
        ]
    output:
        plot = "results/{ICA_path}/sigma_{sigma}/comp_{comp}_gene_length.svg"
    params:
        gene_list_dir = "results/{ICA_path}/sigma_{sigma}/gene_list/"
    conda:
        "../envs/ICA_python.yaml"
    script:
        "../scripts/genomic_overlap/quantifyingUniqueGeneLength.py"


checkpoint defineAdjacentGenes:
    """
        For a gene in a specific annotation, outputs ID of overlapping genes
    """
    input:
        gtf = "data/references/{annotation}.gtf.gz",
        bed = "results/genomic_overlap/gene_bed_merged/{gene}_{annotation}.bed"
    output:
        tsv = "results/genomic_overlap/adjacent_genes/{gene}_{annotation}.tsv"
    conda:
        "../envs/bedtools.yaml"
    script:
        "../scripts/genomic_overlap/defineAdjacentGenes.py"


rule plottingAroundGene:
    """
        From an HGNC ID, plot combined transcripts structure for that gene in
        the three studied annotations, and add all genes that overlap with the
        gene of interest. The gene of interest is defined in green.

        Combined transcripts structure is defined as the shadow of all
        transcripts, without taking into account the coding sequences.
    """
    input:
        beds = partial(getAdjacentGenes, 'bed'),
        HGNC = rna_seq_cartesian_product("data/references/hgnc.txt"),
    output:
        plot = "results/genomic_overlap/plot/{HGNCgene}.svg"
    conda:
        "../envs/ICA_python.yaml"
    script:
        "../scripts/genomic_overlap/plottingAroundGene.py"


rule quantifyRTfromScratch:
    """
        Outputs a three column file containing information about read-throughs
        present in an annotation.

        column1 -> GeneID (read-through)
        column2 -> GeneID (overlapped by read-through)
        column3 -> Number of identical exon between column1 and column2 genes

        To be considered a read-through, the genes needed to overlap at least
        one gene considered in this study, and have at least one identical exon
        with two different genes (and not the same exon for both genes)
    """
    input:
        gtf = "data/references/{annotation}.gtf.gz",
        count_file = "subworkflows/rna_seq_cartesian_product/results/cartesian_product/tissues_NaN.tsv"
    output:
        RT_list = "results/genomic_overlap/readingT/{annotation}.tsv"
    conda:
        "../envs/ICA_python.yaml"
    script:
        "../scripts/genomic_overlap/quantifyRTfromScratch.py"


rule quantifyOverlapfromScratch:
    """
        Outputs a 2-columns tsv containing the gene id and a percentage of exon
        overlap. Sense and antisense are both considered for overlap. Each gene
        is studied inside its own annotation.
    """
    input:
        gtf = "data/references/{annotation}.gtf.gz",
        count_file = "subworkflows/rna_seq_cartesian_product/results/cartesian_product/tissues_NaN.tsv"
    output:
        overlaps = "results/genomic_overlap/overlap/{annotation}.tsv"
    conda:
        "../envs/ICA_python.yaml"
    script:
        "../scripts/genomic_overlap/quantifyOverlapfromScratch.py"


rule plottingAllRT:
    """
        Plots barcharts quantifying the number of genes included in this study
        that are overlapped by read-throughs, separated by genome annotation.
    """
    input:
        RT_list = expand(
            "results/genomic_overlap/readingT/{annotation}.tsv",
            annotation=config["tools"]["annotation"]
        ),
        HGNC = rna_seq_cartesian_product("data/references/hgnc.txt"),
    output:
        plot = "results/genomic_overlap/all_RT.svg"
    conda:
        "../envs/ICA_python.yaml"
    script:
        "../scripts/genomic_overlap/plottingAllRT.py"


rule plottingAllOverlap:
    """
        Violinplots of the percentage of exon overlapped for each genes,
        separated by genome annotation.
    """
    input:
        overlaps = expand(
            "results/genomic_overlap/overlap/{annotation}.tsv",
            annotation=config["tools"]["annotation"]
        ),
        HGNC = rna_seq_cartesian_product("data/references/hgnc.txt")
    output:
        plot = "results/genomic_overlap/all_overlap.svg"
    conda:
        "../envs/ICA_python.yaml"
    script:
        "../scripts/genomic_overlap/plottingAllOverlap.py"


rule getCompBiotypes:
    """
        Outputs a tsv file with biotypes of the significant genes for a comp,
        separated by the side.
    """
    input:
        gtf = "data/references/ensembl98.gtf.gz",
        gene_list = "results/{ICA_path}/sigma_{sigma}/gene_list/comp_{comp}.txt",
        gene_list_tkn = "results/{ICA_path}/sigma_{sigma}/gene_list/.tkn"
    output:
        gene_biotypes = "results/{ICA_path}/sigma_{sigma}/gene_list/comp_{comp}_biotypes.tsv"
    conda:
        "../envs/ICA_python.yaml"
    script:
        "../scripts/genomic_overlap/getCompBiotypes.py"
