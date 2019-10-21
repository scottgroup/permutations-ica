import matplotlib.pyplot as plt
import numpy as np

from gprofiler import GProfiler

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

def plotting(query, fig_path):
    gp = GProfiler(return_dataframe=True)
    df = gp.profile(organism='hsapiens', ordered=True, no_evidences=True, query=query)

    # Keeping only GO terms
    sources = ["GO:BP", "GO:CC", "GO:MF"]
    df = df[df['source'].isin(sources)]

    # Initialisation of the figure
    fig, axes = plt.subplots(nrows=3, ncols=1, sharex=True, figsize=(12, 5))

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
        axes[i].barh(pos, p_val)
        axes[i].set_yticks(pos)
        axes[i].set_yticklabels(label)

        # Formatting plot
        for spine in axes[i].spines.values():
            spine.set_visible(False)

        # Title
        axes[i].set_title(
            source + ' - ' + str(len(source_data)) + ' terms',
            loc='left'
        )

    plt.subplots_adjust(left=0.5, right=0.975, hspace=0.3)
    plt.savefig(fig_path)



# Import genes from file
fpath = snakemake.params.path

# Reading the files
up, down = reading_file(fpath)

# Plotting
plotting(up, snakemake.output.plot_up)
plotting(down, snakemake.output.plot_down)

# glist = [
#     "ENSG00000160882", "ENSG00000089199", "ENSG00000148795", "ENSG00000147465",
#     "ENSG00000203859", "ENSG00000136999", "ENSG00000112936", "ENSG00000085662",
#     "ENSG00000100604", "ENSG00000116133", "ENSG00000214548", "ENSG00000185559",
#     "ENSG00000073060", "ENSG00000138356", "ENSG00000231852", "ENSG00000198886",
#     "ENSG00000137714", "ENSG00000198938", "ENSG00000140459", "ENSG00000210082",
#     "ENSG00000242265", "ENSG00000057252", "ENSG00000023330", "ENSG00000123454",
#     "ENSG00000196205", "ENSG00000198899", "ENSG00000089220", "ENSG00000171951",
#     "ENSG00000143248", "ENSG00000004799", "ENSG00000164125", "ENSG00000143819",
#     "ENSG00000011465", "ENSG00000135111", "ENSG00000106538", "ENSG00000130203",
#     "ENSG00000215187", "ENSG00000123358", "ENSG00000171303", "ENSG00000125730",
#     "ENSG00000144381", "ENSG00000143878", "ENSG00000165449", "ENSG00000138449",
#     "ENSG00000157617", "ENSG00000085563", "ENSG00000198727", "ENSG00000181195",
#     "ENSG00000055732", "ENSG00000075426", "ENSG00000117335", "ENSG00000180176",
#     "ENSG00000179142", "ENSG00000120457", "ENSG00000170962", "ENSG00000072778",
#     "ENSG00000105398", "ENSG00000147256", "ENSG00000100234", "ENSG00000245532",
#     "ENSG00000121966", "ENSG00000134824", "ENSG00000120129", "ENSG00000116741",
#     "ENSG00000145198", "ENSG00000139329", "ENSG00000115461", "ENSG00000188488",
#     "ENSG00000136931", "ENSG00000137801", "ENSG00000142494", "ENSG00000125148",
#     "ENSG00000108654", "ENSG00000111341", "ENSG00000196616", "ENSG00000166148",
#     "ENSG00000111907", "ENSG00000103018", "ENSG00000008394", "ENSG00000162733",
#     "ENSG00000153071", "ENSG00000170899", "ENSG00000072210", "ENSG00000198712",
#     "ENSG00000161513", "ENSG00000236552", "ENSG00000137710", "ENSG00000137463",
#     "ENSG00000133112", "ENSG00000198682", "ENSG00000198840", "ENSG00000089057",
#     "ENSG00000108846", "ENSG00000196352", "ENSG00000023171", "ENSG00000100644",
#     "ENSG00000115828", "ENSG00000063587", "ENSG00000136111", "ENSG00000152661",
#     "ENSG00000162512", "ENSG00000170345", "ENSG00000153234", "ENSG00000119326",
#     "ENSG00000102054", "ENSG00000171724", "ENSG00000164199", "ENSG00000114771",
#     "ENSG00000069667", "ENSG00000130988", "ENSG00000197442", "ENSG00000008283",
#     "ENSG00000198804", "ENSG00000232573", "ENSG00000111269", "ENSG00000082482",
# ]
