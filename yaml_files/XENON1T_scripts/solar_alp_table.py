"""
Make data and tex for solar ALP Bayes factors.
""" 

from log_evidences import *


solar_alp = {}

# Xe column

solar_alp["Xe"] = {"no_3h": bayes_factor(xe1t_alp, xe1t),
                   "3h": bayes_factor(xe1t_3h_alp, xe1t_3h),
                   "bkg_only_3h": bayes_factor(xe1t_alp, xe1t_3h)}

# R-parameter column

no_3h = bayes_factor(alp_r, r)
solar_alp["R"] = {"no_3h": no_3h,
                  "3h": no_3h,
                  "bkg_only_3h": no_3h}

# Xe + R column

solar_alp["Xe + R"] = {"no_3h": bayes_factor(xe1t_alp_r, xe1t_r),
                       "3h": bayes_factor(xe1t_3h_alp_r, xe1t_3h_r),
                       "bkg_only_3h": bayes_factor(xe1t_alp_r, xe1t_3h_r)}

# Xe | R column

solar_alp["Xe | R"] = {k: solar_alp["Xe + R"][k] / solar_alp["R"][k] for k in solar_alp["R"].keys()}


def tex_solar_alp(solar_alp):
    """
    @returns Bayes factors for solar ALP in tex format
    """
    row_names = {"no_3h": "No tritium", "3h": "Tritium", "bkg_only_3h": "Tritium background only"}
    order = ["Xe", "R", "Xe + R", "Xe | R"]
    lines = []

    for k, v in row_names.items():
        line = [v]
        for o in order:
            line.append(r"\num{{{}}}".format(solar_alp[o][k]))

        lines.append(" & ".join(line))

    return "\\\\\n".join(lines)

if __name__ == "__main__":
    print(tex_solar_alp(solar_alp))
