import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
import pandas as pd

from gprofiler import GProfiler

plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams['font.size'] = 6
plt.rcParams['font.family'] = 'monospace'

sources = ["GO:BP", "GO:CC", "GO:MF"]
colors = dict()
colors['down'] = {
    "GO:BP": 'lightcoral',
    "GO:CC": 'r',
    "GO:MF": 'firebrick',
}
colors['up'] = {
    "GO:BP": 'dodgerblue',
    "GO:CC": 'b',
    "GO:MF": 'navy',
}

def get_text(count):
    if count > 1:
        return str(count) + ' genes'
    else:
        return str(count) + ' gene'


def fit_text(text, max_length=37):
    if len(text) > max_length:
        return text[:max_length-3] + '...'
    return text


def mm2inch(*tupl):
    """
    https://stackoverflow.com/questions/14708695/specify-figure-size-in-centimeter-in-matplotlib
    """
    inch = 25.4
    if isinstance(tupl[0], tuple):
        return tuple(i/inch for i in tupl[0])
    else:
        return tuple(i/inch for i in tupl)


def reading_file(path):
    up, down = list(), list()
    with open(path, 'r') as f:
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


def get_GO_df(query):
    gp = GProfiler(return_dataframe=True)
    df = gp.profile(organism='hsapiens', ordered=True, no_evidences=True, query=query)

    # Keeping only GO terms
    df = df[df['source'].isin(sources)]
    df = df[['source', 'native', 'name', 'p_value']]
    df['p_value'] = -np.log10(df['p_value'])
    return df


def plotting_middle():
    weights = sorted(data[comp].values) / max_weigth
    x = np.arange(len(weights))
    lim_down, lim_up = len(down), len(x)-len(up)
    axes[axMid].plot(x[lim_down:lim_up], weights[lim_down:lim_up], 'k')
    axes[axMid].plot(x[:lim_down], weights[:lim_down], 'r')
    axes[axMid].plot(x[lim_up:], weights[lim_up:], 'b')
    axes[axMid].text(x=0.90, y=0.60, s=get_text(len(up)), transform=axes[axMid].transAxes, horizontalalignment='right', verticalalignment='center')
    axes[axMid].text(x=0.10, y=0.40, s=get_text(len(down)), transform=axes[axMid].transAxes, horizontalalignment='left', verticalalignment='center')
    axes[axMid].set_xlabel('Genes')
    axes[axMid].set_ylabel('Weights')
    axes[axMid].axhline(0, color='gray', lw=1, linestyle=':')
    axes[axMid].set_ylim([-1, 1])
    axes[axMid].set_yticks(np.arange(-1, 1.2, step=0.5))
    axes[axMid].yaxis.labelpad = -2
    axes[axMid].xaxis.labelpad = 2


def plotting_GO(df, axes, col, pmin, ori):

    # For every source
    pad = 0.075
    b_pad = 0.02
    r_heigth = (1 - 4*pad - 12*b_pad)/15
    r_width = 0.2

    barsize = 0.004
    pblock_start = 0.75

    for i, source in zip(range(len(sources)), sources):
        source_data = df[df['source'] == source].sort_values('p_value', ascending=False)
        data = source_data[:5]

        # Adding source name
        y_source = 1 - (pad + (pad + 5*r_heigth + 4*b_pad)*i + (4*r_heigth + 2.5*b_pad)/2)
        axes[col].text(x=0, y=y_source, s=source, rotation=90, verticalalignment='center', horizontalalignment='center', weight='bold')
        # Adding decoration
        axes[col].add_patch(patches.Rectangle(
            (0.02, 1 - ((i+1)*(pad + 5*r_heigth + 4*b_pad)-r_heigth)),
            0.01, 5*r_heigth + 4*b_pad,
            color=colors[ori][source]
        ))

        # Adding axis
        axes[col].add_patch(patches.Rectangle(
            (pblock_start, pad), r_width, barsize, color='k', linewidth=barsize/2
        ))
        axes[col].add_patch(patches.Rectangle(
            (pblock_start, pad-0.02), barsize, 0.02, color='k', linewidth=barsize/2
        ))
        axes[col].add_patch(patches.Rectangle(
            (pblock_start+0.20-barsize, pad-0.02), barsize, 0.02, color='k', linewidth=barsize/2
        ))
        axes[col].add_patch(patches.Rectangle(
            (pblock_start+0.10-barsize, pad-0.02), barsize, 0.02, color='k', linewidth=barsize/2
        ))

        axes[col].text(x=pblock_start, y=pad-0.06, s=str(0), horizontalalignment='center')
        axes[col].text(x=pblock_start+0.10-barsize, y=pad-0.06, s=str(int(pmin/2)), horizontalalignment='center')
        axes[col].text(x=pblock_start+0.20-barsize, y=pad-0.06, s=str(pmin), horizontalalignment='center')
        axes[col].text(x=pblock_start+0.2/2, y=pad-0.115, s='-log10(p-value)', horizontalalignment='center')

        for b, data in enumerate(data.values.tolist()):
            source, native, name, p_value = data
            rect_Y = 1 - (pad + i*(pad + 5*r_heigth + 4*b_pad) + b*(r_heigth + b_pad))

            # Writing GO term ID
            char, offset = 37, 0
            if plot_ori == 'up' or plot_ori == 'down':
                char, offset = 64, 0.115
                axes[col].text(
                    x=0.045, y=rect_Y + r_heigth/2, s=fit_text(native),
                    verticalalignment='center'
                )

            # Writing description text
            axes[col].text(
                x=0.045+offset, y=rect_Y + r_heigth/2, s=fit_text(name, char),
                verticalalignment='center'
            )

            # Drawing the p_value rect
            rect_size = r_width/pmin*p_value
            axes[col].add_patch(patches.Rectangle(
                (pblock_start, rect_Y), rect_size, r_heigth,
                color=colors[ori][source]
            ))


# Import genes from file
fpath = snakemake.params.path
comp = snakemake.wildcards.comp

# Reading the files
up, down = reading_file(fpath)

# Reading components_mean from file
data = pd.read_csv(
    snakemake.input.components_mean, sep='\t', index_col=[0,1,2,3]
)
max = np.max(np.max(data))
min = np.min(np.min(data))
max_weigth = np.max([max, -min])


# Initialisation of the figure
plot_ori = snakemake.wildcards.ori
if plot_ori == 'both':
    fig, axes = plt.subplots(
        nrows=1, ncols=4, figsize=mm2inch(178, 65),
        gridspec_kw={'width_ratios':[3,0.1,1,3]}
    )
    axDown, axMid, axUp = 0, 2, 3
    axes[1].axis('off')
elif plot_ori == 'down':
    fig, axes = plt.subplots(
        nrows=1, ncols=3, figsize=mm2inch(178, 65),
        gridspec_kw={'width_ratios':[6,0.1,1]}
    )
    axDown, axMid = 0, 2
    axes[1].axis('off')
elif plot_ori == 'up':
    fig, axes = plt.subplots(
        nrows=1, ncols=3, figsize=mm2inch(178, 65),
        gridspec_kw={'width_ratios':[0.2,1,6]}
    )
    axMid, axUp = 1, 2
    axes[0].axis('off')

# Getting data
down_df = get_GO_df(down)
up_df = get_GO_df(up)
pmin = np.max([np.max(down_df['p_value']), np.max(up_df['p_value'])])
pmin = 90

# Plotting the different parts
plotting_middle()
if plot_ori == 'both' or plot_ori == 'down':
    axes[axDown].axis('off')
    plotting_GO(down_df, axes, axDown, pmin, 'down')

if plot_ori == 'both' or plot_ori == 'up':
    axes[axUp].axis('off')
    plotting_GO(up_df, axes, axUp, pmin, 'up')

fig.suptitle('Expression mode ' + str(int(snakemake.wildcards.comp)+1), fontsize="large", weight='bold')
plt.subplots_adjust(left=0.01, right=0.99, wspace=0.075, top=0.92, bottom=0.125)
plt.savefig(snakemake.output.plot)
