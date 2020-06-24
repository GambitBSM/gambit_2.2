"""
Change prior ranges in yaml.

    python3 make_yaml.py signal.yaml 1e-10 1e-5

gives yaml with prior ranges for the axion couplings between 1e-10 and 1e-5
"""

import yaml
import numpy as np


def alter_couplings(template_name, a, b, suffix=None):

    with open(template_name) as f:
        original = yaml.load(f)

    params = ["gaee", "gan", "gagg"]

    for p in params:
        original["Parameters"]["GeneralALP"][p]["range"] = [float(a), float(b)]

    if suffix is not None:
        original["KeyValues"]["default_output_path"] += suffix

    original["Scanner"]["scanners"]["multinest"]["nlive"] = 1000  # fast repeats

    return original


def alter_tritium(template_name, a, suffix=None):

    with open(template_name) as f:
        original = yaml.load(f)

    original["Parameters"]["XENON1T_NuisanceParameters"]["x_3H"]["sigs"] = [float(a)]
    original["Scanner"]["scanners"]["multinest"]["nlive"] = 1000  # fast repeats

    if suffix is not None:
        original["KeyValues"]["default_output_path"] += suffix

    return original


if __name__ == "__main__":
    
    # Tritium yamls

    sigma = np.linspace(0., 5, 10)

    for i, s in enumerate(sigma):
        d = alter_tritium("bkg_tritium.yaml", s, "_{}".format(i))
        name = "bkg_tritium_{}.yaml".format(i)
        with open(name, "w") as f:
            yaml.dump(d, f, default_flow_style=False)

        d = alter_tritium("signal_tritium.yaml", s, "_{}".format(i))
        name = "signal_tritium_{}.yaml".format(i)
        with open(name, "w") as f:
            yaml.dump(d, f, default_flow_style=False)

    # Signal coupling yamls

    center = -12.5
    width = np.linspace(0., 7.5, 10)

    for i, w in enumerate(width):
        a = center - w
        b = center + w
        d = alter_couplings("signal.yaml", a, b, "_{}".format(i))
        name = "signal_{}.yaml".format(i)
        with open(name, "w") as f:
            yaml.dump(d, f, default_flow_style=False)
