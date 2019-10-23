import matplotlib.pyplot as plt
import numpy as np

from gprofiler import GProfiler

colors = {
    "GO:BP": 'b',
    "GO:CC": 'g',
    "GO:MF": 'r',
}

def reading_file(path):
    up, down = list(), list()
    with open(fpath, 'r') as f:
        for line in f.readlines():
            line = line.strip()
            if line == '>Positive genes':
                reading_up = True
            elif line == '>Negative genes':
                reading_up = False
            else:
                if reading_up:
                    up.append(line)
                else:
                    down.append(line)
    return up, down

def plotting(query, axes, col):
    gp = GProfiler(return_dataframe=True)
    df = gp.profile(organism='hsapiens', ordered=True, no_evidences=True, query=query)

    # Keeping only GO terms
    sources = ["GO:BP", "GO:CC", "GO:MF"]
    df = df[df['source'].isin(sources)]

    # For every source
    for i, source in zip(range(len(sources)), sources):
        source_data = df[df['source'] == source].sort_values('p_value')
        data = source_data[:5]

        # Formatting data
        p_val = list(-np.log(data['p_value']))
        label = data['name'].tolist()
        while len(p_val) < 5:
            p_val.append(0)
            label.append('')

        pos = np.arange(5, 0, -1)

        # Plotting the graph
        axes[i, col].barh(pos, p_val, color=colors[source])
        axes[i, col].set_yticks(pos)
        axes[i, col].set_yticklabels(label)

        # Formatting plot
        for spine in axes[i, col].spines.values():
            spine.set_visible(False)

        # Title
        axes[i, col].set_title(
            source + ' - ' + str(len(source_data)) + ' terms',
            loc='left'
        )

    axes[i, col].set_xlabel('-log(p-value)')



# Import genes from file
fpath = snakemake.params.path

# Reading the files
up, down = reading_file(fpath)

# Initialisation of the figure
fig, axes = plt.subplots(nrows=3, ncols=2, sharex=True, figsize=(15, 4))

# Plotting
plotting(down, axes, 0)
plotting(up, axes, 1)
# plt.subplots_adjust(right=0.975, hspace=0.3)
fig.suptitle('Component ' + snakemake.wildcards.comp, fontsize="x-large")
plt.tight_layout()
plt.savefig(snakemake.output.plot)
