import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from plotting_utils import mm2inch, rcParams
for k, v in rcParams.items():
    plt.rcParams[k] = v

# Loading data
data = pd.read_csv(
    snakemake.input.file, sep='\t',
    index_col=0, header=[0, 1, 2, 3, 4]
)
data = data.groupby(level='aligner', axis=1).mean()

# Plotting
plt.figure(figsize=mm2inch((30,50)))
box = plt.boxplot(
    [
        data['HISAT2'],
        data['STAR'],
        data['tophat2']
    ],
    widths=.8,
    labels=['HISAT2', 'STAR', 'TopHat2']
)

plt.ylim([0,100])
plt.ylabel('% of mate mapped to other chr')
plt.xticks(rotation=90)
plt.tight_layout()
plt.savefig(snakemake.output.plot)
