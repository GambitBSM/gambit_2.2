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
    rc('font', **{'family': 'serif', 'size': 12})
    rc('axes', **{'grid': False, 'titlesize': 12})
    rc('figure', **{'titlesize': 12})
    rc('legend', **{'fontsize': 12, 'title_fontsize': 12, 'handlelength': 1.})
