"""
Styling for plots.
"""

from matplotlib import rc
import matplotlib.pyplot as plt


def make_style():
    """
    Applies style commands to global state
    """
    rc('text', usetex=True)
    rc('text.latex', preamble=r'\usepackage{amsmath}\newcommand{\tritium}{${}^3$H}')
    rc('font', **{'family': 'serif', 'size': 17})
    rc('axes', **{'grid': False, 'titlesize': 15})
    rc('figure', **{'titlesize': 15})
    rc('legend', **{'fontsize': 13, 'title_fontsize': 15, 'handlelength': 1., 'frameon': False})
    rc('xtick', **{'direction': 'in', 'major.size': 7, 'minor.size': 3.5, 'minor.visible': True})
    rc('ytick', **{'direction': 'in', 'major.size': 7, 'minor.size': 3.5, 'minor.visible': True})

def add_logo(fig, x, y, size=0.15, alpha=0.8):
    """
    Adds logo to axis
    """
    im = plt.imread("logo_small.jpg")
    logo_ax = fig.add_axes([x, y, size, size], anchor='NE', zorder=0)
    logo_ax.imshow(im, aspect="equal", alpha=alpha, interpolation="none", zorder=-1)
    logo_ax.axis('off')
