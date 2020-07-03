"""
"""

import numpy as np

# Solar ALP

xe1t_alp = np.array([-21.688189583210804, 2.9673343805847650E-002])
xe1t_alp_r = np.array([-24.378976118255807, 1.8241152697386458E-002])
xe1t_alp_r_wd = np.array([-33.594111494987480, 2.7287332214140699E-002])
xe1t_3h_alp = np.array([-21.414395923422873, 2.7845065229389763E-002])
xe1t_3h_alp_r = np.array([-22.709554198241914, 2.7462782936082987E-002])
xe1t_3h_alp_r_wd = np.array([-32.282489675922214, 3.3312049616034484E-002])
alp_r = np.array([-1.7880299067052574, 1.6627014358437418E-002])
alp_r_wd = np.array([-10.971688775591119, 2.6138961933836878E-002])

# DM ALP

dm_xe1t_alp = np.array([None, None])
dm_xe1t_alp_r_wd = np.array([-32.706501022154001, 4.0278351498150289E-002])
dm_xe1t_3h_alp = np.array([-21.663686916310322, 2.5982589183559394E-002])

# Combine two runs
dm_xe1t_3h_alp_r_wd_1 = np.array([-32.317710884900158, 3.6813797616004028E-002])
dm_xe1t_3h_alp_r_wd_2 = np.array([-32.305534738325065, 3.6688827942681475E-002])
dm_xe1t_3h_alp_r_wd = np.array([0.5 * (dm_xe1t_3h_alp_r_wd_1[0] + dm_xe1t_3h_alp_r_wd_1[0]),
                                (dm_xe1t_3h_alp_r_wd_1[1]**-2 + dm_xe1t_3h_alp_r_wd_1[1]**-2)**-0.5])

# Background only

xe1t = np.array([-22.557779416855389, 6.9891095192451185E-003])
xe1t_3h = np.array([-20.946664875526281, 2.1350060516804027E-002])
# Trivial model - no integration hence no error
r = np.array([-4.7979133853e-01, 0])
r_wd = np.array([-1.1300194552e+01, 0])
# Combine evidences trivially
xe1t_r = np.array([xe1t[0] + r[0], (xe1t[1]**2 + r[1]**2)**0.5])
xe1t_3h_r = np.array([xe1t_3h[0] + r[0], (xe1t_3h[1]**2 + r[1]**2)**0.5])

# Bayes factors

def bayes_factor(r1, r2):
    return np.exp(r1[0] - r2[0])

# Table 3.

t3 = {}

# Table 3. Xe column

t3["Xe"] = {"no_3h": bayes_factor(xe1t_alp, xe1t),
            "3h": bayes_factor(xe1t_3h_alp, xe1t_3h),
            "bkg_only_3h": bayes_factor(xe1t_alp, xe1t_3h)}

# Table 3. R-parameter column

no_3h = bayes_factor(alp_r, r)
t3["R"] = {"no_3h": no_3h,
           "3h": no_3h,
           "bkg_only_3h": no_3h}

# Table 3. Xe + R column

t3["Xe + R"] = {"no_3h": bayes_factor(xe1t_alp_r, xe1t_r),
                "3h": bayes_factor(xe1t_3h_alp_r, xe1t_3h_r),
                "bkg_only_3h": bayes_factor(xe1t_alp_r, xe1t_3h_r)}

# Table 3. Xe | R column

t3["Xe | R"] = {k: t3["Xe + R"][k] / t3["R"][k] for k in t3["R"].keys()}


def tex_t3(t3):
    """
    @returns Bayes factors for table. 3 in tex format
    """
    row_names = {"no_3h": "No tritium", "3h": "Tritium", "bkg_only_3h": "Tritium background only"}
    order = ["Xe", "R", "Xe + R", "Xe | R"]
    lines = []

    for k, v in row_names.items():
        line = [v]
        for o in order:
            line.append(r"\num{{{}}}".format(t3[o][k]))

        lines.append(" & ".join(line))

    return "\\\\\n".join(lines)


if __name__ == "__main__":
    print(tex_t3(t3))
