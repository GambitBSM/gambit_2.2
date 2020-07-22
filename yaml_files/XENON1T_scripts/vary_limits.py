import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm

from style import make_style, add_logo


data_arr = np.loadtxt("xe1t_alp.dat")

width = data_arr[:, 0]
order = width.argsort()
width = width[order]
ln_z_xe1t_alp = data_arr[order, 1]
ln_z_xe1t = -22.557779416855389

bayes_factor = np.exp(ln_z_xe1t_alp - ln_z_xe1t)

make_style()
fig, ax = plt.subplots()
ax_rhs = ax.twinx()

ax.plot(width, bayes_factor, ls="-", marker="o")
ax.set_ylabel("Bayes factor, $B_{10}$")
ax.set_xlabel("Prior width, $w$")
ax.set_title("ALP couplings from $10^{-12.5 - 0.5 w}$ to $10^{-12.5 + 0.5 w}$")

# Limits

ymax = 30.
ax.set_xlim(0, 15)
ax.set_ylim(0, ymax)

# Twin axis showing sigmas

z = [0., 1., 1.5, 1.75]
label = [r"${}\sigma$".format(n) for n in z]
p = norm.sf(z)
y = 1. / p - 1.

ax_rhs.set_ylabel("$Z$-score")
ax_rhs.set_yticks(y)
ax_rhs.set_yticklabels(label)
ax_rhs.set_ylim(0, ymax)
ax_rhs.tick_params(axis='y', which='minor', right=False)

plt.tight_layout()
add_logo(fig, 0.65, 0.725)
plt.savefig("xe1t_alp.pdf")
