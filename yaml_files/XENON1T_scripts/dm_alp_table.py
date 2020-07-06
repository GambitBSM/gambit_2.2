"""
Make data and tex for DM ALP Bayes factors.
""" 

from log_evidences import *


dm_alp = {}

# Xe column

dm_alp["Xe"] = {"no_3h": bayes_factor(dm_xe1t_alp, xe1t),
                "3h": bayes_factor(dm_xe1t_3h_alp, xe1t_3h),
                "bkg_only_3h": bayes_factor(dm_xe1t_alp, xe1t_3h)}

# Xe + R column

dm_alp["Xe + R"] = {"no_3h": bayes_factor(dm_xe1t_alp_r, xe1t_r),
                    "3h": bayes_factor(dm_xe1t_3h_alp_r, xe1t_3h_r),
                    "bkg_only_3h": bayes_factor(dm_xe1t_alp_r, xe1t_3h_r)}
# Xe + R + WD column

dm_alp["Xe + R + WD"] = {"no_3h": bayes_factor(dm_xe1t_alp_r_wd, xe1t_r_wd),
                         "3h": bayes_factor(dm_xe1t_3h_alp_r_wd, xe1t_3h_r_wd),
                         "bkg_only_3h": bayes_factor(dm_xe1t_alp_r_wd, xe1t_3h_r_wd)}

def tex_dm_alp(dm_alp):
    """
    @returns Bayes factors for DM ALP in tex format
    """
    row_names = {"no_3h": "No tritium", "3h": "Tritium", "bkg_only_3h": "Tritium background only"}
    order = ["Xe", "Xe + R", "Xe + R + WD"]
    lines = []

    for k, v in row_names.items():
        line = [v]
        for o in order:
            line.append(r"\num{{{}}}".format(dm_alp[o][k]))

        lines.append(" & ".join(line))

    return "\\\\\n".join(lines)

if __name__ == "__main__":
    print(tex_dm_alp(dm_alp))
