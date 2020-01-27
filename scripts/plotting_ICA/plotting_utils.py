
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
