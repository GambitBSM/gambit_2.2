import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm

from style import make_style


data_arr = np.loadtxt("xe1t_3h.dat")
width = data_arr[:, 0]
order = width.argsort()
width = width[order]
ln_z_xe1t_3h = data_arr[order, 1]

data_arr = np.loadtxt("xe1t_3h_alp.dat")
width = data_arr[:, 0]
order = width.argsort()
width = width[order]
ln_z_xe1t_3h_alp = data_arr[order, 1]

ln_z_xe1t = -22.557779416855389

alp_3h_vs_bkg_3h = np.exp(ln_z_xe1t_3h_alp - ln_z_xe1t_3h)
bkg_3h_vs_bkg = np.exp(ln_z_xe1t_3h - ln_z_xe1t)

make_style()
fig, ax = plt.subplots()
ax_rhs = ax.twinx()
col_rhs = "red"

ax.plot(width, alp_3h_vs_bkg_3h, ls="-", marker="o", label="Signal + tritium versus background + tritium")
ax.plot(width, bkg_3h_vs_bkg, ls="-", c="tab:green", marker="o", label="Background + tritium versus background")
ax.legend()
ax.set_title(r"Tritium component $\log_{10}\frac{\alpha_t}{\text{mol/mol}} = -27 \pm \sigma$")
ax.set_ylabel(r"Bayes factor, $B$")
ax.set_xlabel(r"Uncertainty in tritium component, $\sigma$")

# Limits

ymax = 7.

ax.set_xlim(0, 5.5)
ax.set_ylim(0, ymax)

# Twin axis showing sigmas

z = [0., 0.5, 1., 1.5, 1.75]
label = [r"${}\sigma$".format(n) for n in z]
p = norm.sf(z)
y = 1. / p - 1.

ax_rhs.set_ylabel("$Z$-score", color=col_rhs)
ax_rhs.spines['right'].set_color(col_rhs)
ax_rhs.tick_params(axis='x', colors=col_rhs)
ax_rhs.tick_params(axis='y', colors=col_rhs)
ax_rhs.set_yticks(y)
ax_rhs.set_yticklabels(label)
ax_rhs.set_ylim(0, ymax)

plt.tight_layout()
plt.savefig("xe1t_3h_alp.pdf")
