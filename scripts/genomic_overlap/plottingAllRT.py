import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from plotting_utils import mm2inch, rcParams
for k, v in rcParams.items():
    plt.rcParams[k] = v


n_RT = list()
for overRT_path in snakemake.input.RT_list:
    annotation = overRT_path.split('/')[-1].split('.')[0]
    overRT = pd.read_csv(overRT_path, sep='\t')

    n_overlap = overRT.groupby('gene').count().reset_index()
    read_through = n_overlap[n_overlap['overlap'] > 1]
    n_RT.append((annotation, len(read_through)))

data = pd.DataFrame(n_RT, columns=['variable', 'value'])

# Plotting violin
fig = plt.figure(figsize=mm2inch((80, 45)))
sns.barplot(y='variable', x='value', data=data)
plt.xlabel('N genes overlapped by read-through')
plt.ylabel('')
plt.tight_layout()
plt.savefig(snakemake.output.plot)
