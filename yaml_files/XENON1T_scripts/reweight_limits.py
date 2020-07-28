"""
Re-weight the evidence to adjust breadth of ALP couplings prior
===============================================================
"""

import numpy as np
import log_evidences_arxiv_v2 as ln_z


# original prior range
ya = -20.
yb = -3.


def unit_breadth(a, b):
    """
    @returns Cuts on unit interval
    """
    ua = (a - ya) / (yb - ya)
    ub = (b - ya) / (yb - ya)
    return (ua, ub)

def norm_factor(a, b):
    """
    @returns Factor that adjuss normalization of prior
    """
    return ((yb - ya) / (b - a))**3


if __name__ == "__main__":

    name = "/home/andrew/Desktop/v2/xe1t_alp/scanner_plugins/MultiNest/native-.txt"
    data = np.loadtxt(name)

    center = -12.5
    widths = np.linspace(0., 15., 15)

    for width in widths:
        a = center - 0.5 * width
        b = center + 0.5 * width
        ua, ub = unit_breadth(a, b)
        where = np.logical_and.reduce((data[:, 2] < ub, data[:, 2] > ua,
                                       data[:, 3] < ub, data[:, 3] > ua,
                                       data[:, 4] < ub, data[:, 4] > ua))
        weight = data[where, 0].sum()
        log_evidence = ln_z.xe1t_alp[0] + np.log(norm_factor(a, b)) + np.log(weight)
        print(width, log_evidence)
