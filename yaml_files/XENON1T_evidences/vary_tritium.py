import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm

from style import make_style


data_arr_bkg = np.loadtxt("vary_tritium_bkg.dat")
width = data_arr[:, 0]
order = width.argsort()
width = width[order]
log_z_bkg = data_arr[order, 1]

data_arr_signal = np.loadtxt("vary_tritium_signal.dat")
width = data_arr[:, 0]
order = width.argsort()
width = width[order]
log_z_signal = data_arr[order, 1]

bayes_factor = np.exp(log_z_signal - log_z_bkg)

make_style()
fig, ax = plt.subplots()
ax_rhs = ax.twinx()
col_rhs = "red"

ax.plot(width, bayes_factor, ls="-", marker="o")
plt.title(r"Prior $\log_{10} \alpha_t = -27 \pm \sigma$")
ax.set_ylabel("Bayes factor, $B_{10}$")
ax.set_xlabel("$\sigma$")
ax.set_xlim(0, 15)
ax.set_ylim(0, 30)

z = [0., 0.5, 1., 1.5, 1.75]
p = norm.sf(z)
y = 1. / p - 1.

ax_rhs.set_ylabel("Bayesian significance, $Z$", color=col_rhs)
ax_rhs.spines['right'].set_color(col_rhs)
ax_rhs.tick_params(axis='x', colors=col_rhs)
ax_rhs.tick_params(axis='y', colors=col_rhs)
ax_rhs.set_yticklabels(z)
ax_rhs.set_yticks(y)
ax_rhs.set_ylim(0, 30)

plt.savefig("vary_tritium.pdf")
