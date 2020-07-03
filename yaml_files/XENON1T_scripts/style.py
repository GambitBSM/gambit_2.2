"""
Styling for plots.
"""

from matplotlib import rc

def make_style():
    """
    Applies style commands to global state
    """
    rc('text', usetex=True)
    rc('text.latex', preamble=r'\usepackage{amsmath}')
    rc('font', **{'family': 'serif', 'size': 14})
    rc('axes', **{'grid': False, 'titlesize': 14})
    rc('figure', **{'titlesize': 14})
    rc('legend', **{'fontsize': 14, 'title_fontsize': 14, 'handlelength': 1.})
