"""
Make data and tex for DM ALP Bayes factors.
"""

import log_evidences_arxiv_v2 as logz
from bayes import bayes_factor


dm_alp = {}

# Xe column

dm_alp["Xe"] = {"no_3h": bayes_factor(logz.dm_xe1t_alp, logz.xe1t),
                "3h": bayes_factor(logz.dm_xe1t_3h_alp, logz.xe1t_3h),
                "bkg_only_3h": bayes_factor(logz.dm_xe1t_alp, logz.xe1t_3h)}

# R-parameter column

no_3h = bayes_factor(logz.dm_alp_r, logz.r)
dm_alp["R"] = {"no_3h": no_3h,
               "3h": no_3h,
               "bkg_only_3h": no_3h}

# Xe + R column

dm_alp["Xe + R"] = {"no_3h": bayes_factor(logz.dm_xe1t_alp_r, logz.xe1t_r),
                    "3h": bayes_factor(logz.dm_xe1t_3h_alp_r, logz.xe1t_3h_r),
                    "bkg_only_3h": bayes_factor(logz.dm_xe1t_alp_r, logz.xe1t_3h_r)}
# Xe + R + WD column

dm_alp["Xe + R + WD"] = {"no_3h": bayes_factor(logz.dm_xe1t_alp_r_wd, logz.xe1t_r_wd),
                         "3h": bayes_factor(logz.dm_xe1t_3h_alp_r_wd, logz.xe1t_3h_r_wd),
                         "bkg_only_3h": bayes_factor(logz.dm_xe1t_alp_r_wd, logz.xe1t_3h_r_wd)}

# Xe | R column

dm_alp["Xe | R"] = {k: dm_alp["Xe + R"][k] / dm_alp["R"][k] for k in dm_alp["R"].keys()}

# R + WD column

dm_alp["R + WD"] = {"no_3h": bayes_factor(logz.dm_alp_r_wd, logz.r_wd),
                    "3h": bayes_factor(logz.dm_alp_r_wd, logz.r_wd),
                    "bkg_only_3h": bayes_factor(logz.dm_alp_r_wd, logz.r_wd)}

# Xe | R + WD column

dm_alp["Xe | R + WD"] = {k: dm_alp["Xe + R + WD"][k] / dm_alp["R + WD"][k] for k in dm_alp["R"].keys()}

def tex_dm_alp(dm_alp):
    """
    @returns Bayes factors for DM ALP in tex format
    """
    row_names = {"no_3h": "No \ce{^3H}", "3h": "\ce{^3H}", "bkg_only_3h": "\ce{^3H} background only"}
    order = ["Xe", "Xe + R", "Xe + R + WD", "Xe | R", "Xe | R + WD"]
    lines = []

    for k, v in row_names.items():
        line = [v]
        for o in order:
            line.append(r"\num{{{}}}".format(dm_alp[o][k]))

        lines.append(" & ".join(line))

    return "\\\\\n".join(lines)

if __name__ == "__main__":
    print(tex_dm_alp(dm_alp))

    # Check dependence on eta prior
    flat_eta = bayes_factor(logz.dm_xe1t_alp_r_wd, logz.xe1t_r_wd)
    log_eta = bayes_factor(logz.dm_xe1t_alp_r_wd_logeta, logz.xe1t_r_wd)
    print("Bayes factor with log eta prior", log_eta)
    print("corresponding to an increase by factor", log_eta / flat_eta)
