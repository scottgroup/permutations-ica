import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Loading data
data = pd.read_csv(
    snakemake.input.data, names=['IDs', 'metric'],
    sep='\t', header=0
)

plt.plot(data['IDs'], data['metric'])
plt.show()
