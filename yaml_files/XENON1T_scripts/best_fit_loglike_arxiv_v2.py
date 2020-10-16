"""
Best-fit loglikes for version 2 on arXiv

Compared to version 1, we corrected a bug in the inverse Primakoff contributions.
"""

from best_fit_loglike_arxiv_v1 import xe1t, xe1t_r, xe1t_r_wd, xe1t_3h, xe1t_3h_r, xe1t_3h_r_wd

# Solar ALP
xe1t_alp = -14.63018597625
xe1t_alp_r = -21.7191121066993
xe1t_alp_r_wd = -27.7529764920372

# Solar ALP with tritium
xe1t_3h_alp = -14.63293792966
xe1t_3h_alp_r = -16.81193382537
xe1t_3h_alp_r_wd = -22.85894817789

# DM ALP
dm_xe1t_alp = -13.6098820469
dm_xe1t_alp_r = -13.6271472324
dm_xe1t_alp_r_wd = -21.7292624563

# DM ALP with tritium
dm_xe1t_3h_alp = -12.9739359185
dm_xe1t_3h_alp_r = -13.0127377906
dm_xe1t_3h_alp_r_wd = -21.0908853881
