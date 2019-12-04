import os

import numpy as np
import pandas as pd
from pathlib import Path
import pysam
import pybedtools

from get_gene_gtf import get_gene
from get_gene_gtf import get_bam_variables


# Importing data
bams = snakemake.input.bams
df_gene, start, stop, strand, chr = get_gene(
    snakemake.input.gtf, snakemake.wildcards.gene
)

# Importing datasets depth
datasets_depth = pd.read_csv(
    snakemake.input.datasets_depth,
    sep='\t', index_col=[0]
)['size'].to_dict()

# Creating count matrix and header list
bed_mat = np.zeros((stop-start, len(bams)))
header = list()

for i, bam in enumerate(bams):
    # Get variables
    bam_var = get_bam_variables(bam)
    header.append(bam_var)

    # Load BAM file
    samfile = pysam.AlignmentFile(bam, "rb")

    # Creating gene BAM file
    basename = snakemake.wildcards.gene + '.bam'
    temp_fname = os.path.join(os.path.dirname(bam), basename)
    with pysam.Samfile(temp_fname, 'wb', template=samfile) as fout:
        for read in samfile.fetch(str(chr), start, stop):
            fout.write(read)
    samfile.close()

    # Reading tempBam file with pyBedTools
    bed = pybedtools.BedTool(temp_fname)
    bed_cov = bed.genome_coverage(
        bga=True, split=True, strand=strand
    )

    for x in bed_cov:
        _start = max(x.start, start)
        _stop = min(x.stop, stop)

        if _start <= stop and _stop > start:
            _start, _stop = _start-start, _stop-start
            bed_mat[_start:_stop, i] += x.count

    # Add scaling factor for dataset depth
    bed_mat[:, i] *= float(datasets_depth[bam_var['dataset']])

# Touching the gene token file
Path(snakemake.output.gene_BAM_tkn).touch()

# Create output dataframe
df_header = pd.DataFrame(header)
df = pd.DataFrame(bed_mat, columns=pd.MultiIndex.from_frame(df_header))
df.to_csv(snakemake.output.gene_coverage, sep='\t')
