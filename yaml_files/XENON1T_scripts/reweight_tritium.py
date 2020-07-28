"""
Re-weight the evidence to adjust sigma in alpha_t prior
========================================================
"""

import numpy as np
import log_evidences_arxiv_v2 as ln_z
from scipy.stats import norm


# original prior
sigma = 3.
mu = -27.


def reweight_factors(u, width):
    """
    @returns Re-weight factor
    """
    log_alpha = norm(-27., 3.).ppf(u)
    z = ((log_alpha - mu) / sigma)**2
    new = (sigma / width)**2 * z
    return sigma / width * np.exp(-0.5 * (new - z))
    

if __name__ == "__main__":

    widths = np.linspace(0.5, 5, 20)
    name = "/home/andrew/Desktop/v2/xe1t_3h_alp/scanner_plugins/MultiNest/native-.txt"
    col = 7  # 4 for background and 7 for solar ALP
    log_evidence_original = ln_z.xe1t_3h_alp[0]  # ln_z.xe1t_3h_alp[0] for solar ALP
    data = np.loadtxt(name)
    weights = data[:, 0]
    u = data[:, col]

    for width in widths:
        reweight_factors_ = reweight_factors(u, width)
        weight_sum = (reweight_factors_ * weights).sum()
        log_evidence = log_evidence_original + np.log(weight_sum)
        print(width, log_evidence)
