import os
import pandas as pd

from pathlib import Path

# Importing data
filtered_genes = pd.read_csv(
    snakemake.input.filt_genes, sep='\t', index_col=[0, 1, 2, 3]
)

# For every component
for comp in filtered_genes.columns:
    _genes = filtered_genes[comp]

    # Selecting and ranking genes
    _genes_up = _genes[_genes > 0].sort_values(ascending=False)
    _genes_list_up = _genes_up.index.get_level_values('ensembl_id').tolist()

    _genes_down = _genes[_genes < 0].sort_values(ascending=True)
    _genes_list_down = _genes_down.index.get_level_values('ensembl_id').tolist()

    # Writing genes to file
    os.makedirs(snakemake.params.directory, exist_ok=True)
    file_path = os.path.join(
        snakemake.params.directory,
        "comp_{comp}.txt".format(comp=comp)
    )

    with open(file_path, 'w') as f:
        f.write('>Positive genes\n')
        for gene in _genes_list_up:
            f.write(gene + '\n')

        f.write('>Negative genes\n')
        for gene in _genes_list_down:
            f.write(gene + '\n')

Path(snakemake.output.tkn).touch(exist_ok=True)
