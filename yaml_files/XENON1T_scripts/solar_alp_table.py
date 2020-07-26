"""
Make data and tex for solar ALP Bayes factors.
"""

import log_evidences_arxiv_v2 as logz
from bayes import bayes_factor, z_score


solar_alp = {}

# Xe column

solar_alp["Xe"] = {"no_3h": bayes_factor(logz.xe1t_alp, logz.xe1t),
                   "3h": bayes_factor(logz.xe1t_3h_alp, logz.xe1t_3h),
                   "bkg_only_3h": bayes_factor(logz.xe1t_alp, logz.xe1t_3h)}

# R-parameter column

no_3h = bayes_factor(logz.alp_r, logz.r)
solar_alp["R"] = {"no_3h": no_3h,
                  "3h": no_3h,
                  "bkg_only_3h": no_3h}

# Xe + R column

solar_alp["Xe + R"] = {"no_3h": bayes_factor(logz.xe1t_alp_r, logz.xe1t_r),
                       "3h": bayes_factor(logz.xe1t_3h_alp_r, logz.xe1t_3h_r),
                       "bkg_only_3h": bayes_factor(logz.xe1t_alp_r, logz.xe1t_3h_r)}
# Xe + R + WD column

solar_alp["Xe + R + WD"] = {"no_3h": bayes_factor(logz.xe1t_alp_r_wd, logz.xe1t_r_wd),
                            "3h": bayes_factor(logz.xe1t_3h_alp_r_wd, logz.xe1t_3h_r_wd),
                            "bkg_only_3h": bayes_factor(logz.xe1t_alp_r_wd, logz.xe1t_3h_r_wd)}

# Xe | R column

solar_alp["Xe | R"] = {k: solar_alp["Xe + R"][k] / solar_alp["R"][k] for k in solar_alp["R"].keys()}

# R + WD column

solar_alp["R + WD"] = {"no_3h": bayes_factor(logz.alp_r_wd, logz.r_wd),
                       "3h": bayes_factor(logz.alp_r_wd, logz.r_wd),
                       "bkg_only_3h": bayes_factor(logz.alp_r_wd, logz.r_wd)}

# Xe | R + WD column

solar_alp["Xe | R + WD"] = {k: solar_alp["Xe + R + WD"][k] / solar_alp["R + WD"][k] for k in solar_alp["R"].keys()}


def tex_solar_alp(solar_alp):
    """
    @returns Bayes factors for solar ALP in tex format
    """
    row_names = {"no_3h": "No \ce{^3H}", "3h": "\ce{^3H}", "bkg_only_3h": "\ce{^3H} background only"}
    order = ["Xe", "Xe + R", "Xe + R + WD", "Xe | R", "Xe | R + WD"]
    lines = []

    for k, v in row_names.items():
        line = [v]
        for o in order:
            line.append(r"\num{{{}}}".format(solar_alp[o][k]))

        lines.append(" & ".join(line))

    return "\\\\\n".join(lines)

if __name__ == "__main__":
    print(tex_solar_alp(solar_alp))
    print("z_score_for_text", z_score(logz.xe1t_alp, logz.xe1t))
    print("R + WD", bayes_factor(logz.alp_r_wd, logz.r_wd))
    print("xenon + R + WD", bayes_factor(logz.xe1t_alp_r_wd, logz.xe1t_r_wd))
    print("xenon | R + WD", bayes_factor(logz.xe1t_alp_r_wd, logz.xe1t_r_wd) / bayes_factor(logz.alp_r_wd, logz.r_wd))
