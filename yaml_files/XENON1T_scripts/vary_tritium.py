import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm

from style import make_style, add_logo


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

ax.plot(width, alp_3h_vs_bkg_3h, ls="-", marker="o", label=r"Solar ALP + \tritium~vs.~background + \tritium")
ax.plot(width, bkg_3h_vs_bkg, ls="-", c="tab:green", marker="o", label=r"Background + \tritium~vs.~background")
ax.legend(loc="upper left")
ax.set_title(r"\tritium~component $\log_{10}\frac{\alpha_t}{\text{mol/mol}} = -27 \pm \sigma$")
ax.set_ylabel(r"Bayes factor, $B$")
ax.set_xlabel(r"Uncertainty in \tritium~component, $\sigma$")

# Limits

ymax = 7.

ax.set_xlim(0, 5.5)
ax.set_ylim(0, ymax)

# Twin axis showing sigmas

z = [0., 0.5, 1., 1.5, 1.75]
label = [r"${}\sigma$".format(n) for n in z]
p = norm.sf(z)
y = 1. / p - 1.

ax_rhs.set_ylabel("$Z$-score")
ax_rhs.set_yticks(y)
ax_rhs.set_yticklabels(label)
ax_rhs.set_ylim(0, ymax)
ax_rhs.tick_params(axis='y', which='minor', right=False)

plt.tight_layout()
# add_logo(fig, 0.67, 0.71)
plt.savefig("xe1t_3h_alp.pdf")
