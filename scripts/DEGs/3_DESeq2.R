# Script from : https://github.com/snakemake-workflows/rna-seq-star-deseq2/blob/master/scripts/deseq2.R
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

dds <- readRDS(snakemake@input[["rds"]])

# Defining wildcards
variable <- snakemake@wildcards[["variable"]]
tool1 <- snakemake@wildcards[["tool"]]
tool2 <- snakemake@wildcards[["tool2"]]

# Calculating contrast
res <- results(dds, contrast=c(variable,tool1,tool2), parallel=parallel)

# Writing to file
write.csv(
    as.data.frame(res),
    file=snakemake@output[["results"]],
    quote=FALSE
)
