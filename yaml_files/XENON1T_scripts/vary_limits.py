"""
Two-dimensional plot showing Bayes factor sensitivity to prior center and width
===============================================================================
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors

from style import make_style, add_logo
import log_evidences_arxiv_v2 as ln_z


data_arr = np.loadtxt("xe1t_alp.dat")
center = data_arr[:, 0]
width = data_arr[:, 1]
ln_z_xe1t_alp = data_arr[:, 2]
bayes_factor = np.exp(ln_z_xe1t_alp - ln_z.xe1t[0])

make_style()
fig, ax = plt.subplots()

x, y = np.meshgrid(np.unique(width), np.unique(center))
z = bayes_factor.reshape(len(x), len(y))
z = np.ma.masked_invalid(z)

log_levels = np.arange(-1., 3.1, 0.125)
levels = np.power(10., log_levels)
cf = ax.contourf(x, y, z.T, levels, norm=colors.LogNorm(), extend="both", cmap='hot_r')

cbar = fig.colorbar(cf)
cbar.set_ticks([1e-1, 1, 1e1, 1e2, 1e3])
cbar.set_label("Bayes factor, $B_{10}$")
ax.set_xlabel("Prior width, $w$")
ax.set_ylabel("Prior center, $c$")
ax.set_title("Solar ALP couplings from $10^{c - 0.5 w}$ to $10^{c + 0.5 w}$")

p = ax.scatter(17., -11.5, color="royalblue", s=50, marker="H", zorder=5, label="Priors used in analysis")
p.set_clip_on(False)

ax.legend()

# Limits

#ax.set_ylim(-13.5, -9.5)
ax.set_xlim(0., 17.)

plt.tight_layout()
# add_logo(fig, 0.65, 0.725)
plt.savefig("xe1t_alp.pdf")
