# Lines from https://github.com/snakemake-workflows/rna-seq-star-deseq2/
log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library("DESeq2")

parallel <- FALSE
if (snakemake@threads > 1) {
    library("BiocParallel")
    # setup parallelization
    register(MulticoreParam(snakemake@threads))
    parallel <- TRUE
}

# Loading count matrix from file, converting to integer
counts <- as.matrix(
    read.table(
        snakemake@input[["counts"]], header=TRUE,
        row.names="symbol", check.names=FALSE
    )
)
mode(counts) <- "integer"

# Loading samples information from file.
samples <- read.table(
    snakemake@input[["samples"]], header=TRUE,
    row.names="sample", check.names=FALSE
)

# Calculating DESeq2
dds <- DESeqDataSetFromMatrix(
    countData=counts,
    colData=samples,
    design= as.formula(paste("~", snakemake@wildcards[["variable"]]))
)

# Defining the base level
dds[[snakemake@wildcards[["variable"]]]] <- relevel(
    dds[[snakemake@wildcards[["variable"]]]],
    ref=snakemake@wildcards[["tool"]]
)

dds <- DESeq(dds, parallel=parallel)

saveRDS(dds, file=snakemake@output[["rds"]])
