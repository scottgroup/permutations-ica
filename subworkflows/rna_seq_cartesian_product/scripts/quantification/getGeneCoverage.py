import os

import numpy as np
import pandas as pd
import pysam
import pybedtools

from get_gene_gtf import get_gene


def get_bam_variables(fname):
    """ """
    bam_var = dict()

    for tool, vars in snakemake.config['tools'].items():
        for var in vars:
            if var in fname.split('/'):
                bam_var[tool] = var

    for tissue, datasets in snakemake.config['datasets'].items():
        if tissue in fname:
            bam_var['tissue'] = tissue
            for dataset in datasets:
                if dataset in fname:
                    bam_var['dataset'] = dataset

    return bam_var

# Importing data
params = eval(snakemake.params.__str__())
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

    # Creating temporary bam file
    basename = str(start) + '_' + str(stop) + '_' + strand + '.bam'
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

    # Removing temp file
    os.remove(temp_fname)

# Create output dataframe
df_header = pd.DataFrame(header)
df = pd.DataFrame(bed_mat, columns=pd.MultiIndex.from_frame(df_header))
df.to_csv(snakemake.output.gene_coverage, sep='\t')
