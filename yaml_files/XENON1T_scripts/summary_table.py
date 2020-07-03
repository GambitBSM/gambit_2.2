"""
Make data and tex for summary of DM and solar ALP models.
""" 

from log_evidences import *


summary = {}

# Xe column

summary["Xe"] = {"no_3h": evidence(xe1t),
                 "3h": evidence(xe1t_3h),
                 "alp": evidence(xe1t_alp),
                 "alp_3h": evidence(xe1t_3h_alp),
                 "dm_alp": evidence(dm_xe1t_alp),
                 "dm_alp_3h": evidence(dm_xe1t_3h_alp)}

# Xe + R column

summary["Xe + R"] = {"no_3h": evidence(xe1t_r),
                     "3h": evidence(xe1t_3h_r),
                     "alp": evidence(xe1t_alp_r),
                     "alp_3h": evidence(xe1t_3h_alp_r),
                     "dm_alp": evidence(dm_xe1t_alp_r),
                     "dm_alp_3h": evidence(dm_xe1t_3h_alp_r)}

# Xe + R + WD column

summary["Xe + R + WD"] = {"no_3h": evidence(xe1t_r_wd),
                          "3h": evidence(xe1t_3h_r_wd),
                          "alp": evidence(xe1t_alp_r_wd),
                          "alp_3h": evidence(xe1t_3h_alp_r_wd),
                          "dm_alp": evidence(dm_xe1t_alp_r_wd),
                          "dm_alp_3h": evidence(dm_xe1t_3h_alp_r_wd)}

# Greatest evidences for each data set

best = {k: max(v.values()) for k, v in summary.items()}

def tex_summary(summary):
    """
    @returns Summary in tex format
    """
    row_names = {"no_3h": "Background", "3h": "Background + \ce{^3H}",
                 "alp": "Solar ALP", "alp_3h": "Solar ALP + \ce{^3H}",
                 "dm_alp": "DM ALP", "dm_alp_3h": "DM ALP + \ce{^3H}"}
    order = ["Xe", "Xe + R", "Xe + R + WD"]
    lines = []

    for k, v in row_names.items():
        line = [v]
        for o in order:
            line.append(r"\num{{{}}}".format(summary[o][k] / best[o]))

        lines.append(" & ".join(line))

    return "\\\\\n".join(lines)

if __name__ == "__main__":
    print(tex_summary(summary))
