"""
Change prior ranges in yaml.

    python3 make_yaml.py signal.yaml 1e-10 1e-5

gives yaml with prior ranges for the axion couplings between 1e-10 and 1e-5
"""

import yaml


def make_yaml(template_name, a, b):

    with open(template_name) as f:
        original = yaml.load(f)

    params = ["gaee", "gan", "gagg"]

    for p in params:
        original["Parameters"]["GeneralALP"][p]["range"] = [a, b]

    return yaml.dump(original, default_flow_style=False)

    
if __name__ == "__main__":
    
    import sys

    template_name = sys.argv[1]
    a = float(sys.argv[2])
    b = float(sys.argv[3])

    print(make_yaml(template_name, a, b))
