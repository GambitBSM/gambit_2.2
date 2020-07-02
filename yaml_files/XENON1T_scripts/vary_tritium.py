import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm

from style import make_style


data_arr_bkg = np.loadtxt("vary_tritium_bkg.dat")
width = data_arr[:, 0]
order = width.argsort()
width = width[order]
log_z_bkg_tritium = data_arr[order, 1]

data_arr_signal = np.loadtxt("vary_tritium_signal.dat")
width = data_arr[:, 0]
order = width.argsort()
width = width[order]
log_z_signal_tritium = data_arr[order, 1]

log_z_bkg = -22.761451395014724

bf_signal_t_vs_bkg_t = np.exp(log_z_signal_tritium - log_z_bkg_tritium)
bf_bkg_t_vs_bkg = np.exp(log_z_bkg_tritium - log_z_bkg)

make_style()
fig, ax = plt.subplots()
ax_rhs = ax.twinx()
col_rhs = "red"

ax.plot(width, bf_signal_t_vs_bkg_t, ls="-", marker="o", label="Signal + tritium versus background + tritium")
ax.plot(width, bf_bkg_t_vs_bkg, ls="-", marker="o", label="Background + tritium versus background")
ax.legend()
plt.title(r"Prior $\log_{10} \alpha_t = -27 \pm \sigma$")
ax.set_ylabel("Bayes factor, $B$")
ax.set_xlabel("Uncertainty in tritium component, $\sigma$")

# Limits

ymax = 30.

ax.set_xlim(0, 15)
ax.set_ylim(0, ymax)

# Twin axis showing sigmas

z = [0., 0.5, 1., 1.5, 1.75]
p = norm.sf(z)
y = 1. / p - 1.

ax_rhs.set_ylabel("Bayesian significance, $Z$", color=col_rhs)
ax_rhs.spines['right'].set_color(col_rhs)
ax_rhs.tick_params(axis='x', colors=col_rhs)
ax_rhs.tick_params(axis='y', colors=col_rhs)
ax_rhs.set_yticklabels(z)
ax_rhs.set_yticks(y)
ax_rhs.set_ylim(0, ymax)

plt.savefig("vary_tritium.pdf")
