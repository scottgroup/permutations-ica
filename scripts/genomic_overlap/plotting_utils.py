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


def mm2inch(*tupl):
    """
    https://stackoverflow.com/questions/14708695/specify-figure-size-in-centimeter-in-matplotlib
    """
    inch = 25.4
    if isinstance(tupl[0], tuple):
        return tuple(i/inch for i in tupl[0])
    else:
        return tuple(i/inch for i in tupl)

# Figure params
rcParams = {
    'svg.fonttype': 'none',
    'font.size': 6,
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial']
}

# Aligner colors
colors = {
    'tophat2': '#2ca02cff',
    'HISAT2': '#1f77b4ff',
    'STAR': '#ff7f0eff',
    'refseq': '#2ca02cff',
    'ensembl92': '#1f77b4ff',
    'ensembl98': '#ff7f0eff'
}
