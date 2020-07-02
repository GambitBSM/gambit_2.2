"""
Change prior ranges in yaml.

    python3 make_yaml.py

gives you lots of yamls I want to run.
"""

import yaml
import numpy as np
import sys
import os


def alter_couplings(template_name, a, b, suffix=None):

    with open(template_name) as f:
        original = yaml.load(f)

    params = ["gaee", "gaN", "gagg"]

    for p in params:
        original["Parameters"]["GeneralALP"][p]["range"] = [float(a), float(b)]

    if suffix is not None:
        original["KeyValues"]["default_output_path"] += suffix

    # don't adjust nlive - otherwise get noisy result
    original["Scanner"]["scanners"]["multinest"]["efr"] = 0.8  # fast repeats

    return original


def alter_tritium(template_name, a, suffix=None):

    with open(template_name) as f:
        original = yaml.load(f)

    original["Parameters"]["XENON1T_NuisanceParameters"]["x_3H"]["sigs"] = [float(a)]
    original["Scanner"]["scanners"]["multinest"]["nlive"] = 1000  # fast repeats
    original["Scanner"]["scanners"]["multinest"]["efr"] = 0.8  # fast repeats

    if suffix is not None:
        original["KeyValues"]["default_output_path"] += suffix

    return original


if __name__ == "__main__":

    if len(sys.argv) < 2:
        print "Error! Missing parent file.\nUsage make_yaml.py <parent_file>\n"
        exit()

    # Parent file
    filename = sys.argv[1]

    # What I am varying, width in ALP model or tritium sigma on XENON1T nuisance model
    alp = False
    tritium = False
    if "alp" in filename and not "3h" in filename:
      alp = True
    elif "3h" in filename:
      tritium = True

    # Pattern for daughter files
    daughter = '.'.join(filename.split('.')[:-1])
    extension = '.' + filename.split('.')[-1]

    # Make all daughters in a separate directory
    if not os.path.isdir(daughter):
      os.mkdir(daughter)

    # Create temp file with all !imports resolved
    tempfile = filename + "_temp"
    with open(filename) as fin, open(tempfile, 'w') as fout:
        for line in fin:

          if "!import" in line: # if there is a !include, resolve it
             includefile = line.split('import')[-1].strip()
             fout.write(line.split('!import')[:-1][0] + '\n')
             with open(includefile) as finclude:
                 for line2 in finclude:
                     fout.write(line2)
         
          else: # otherwise just copy it
              fout.write(line)

    if tritium:

        # Tritium yamls

        sigma = np.linspace(0.5, 5, 20)

        for i, s in enumerate(sigma):
            suffix = "_{}".format(i)
            d = alter_tritium(tempfile, s, suffix)

            name = daughter + '/' + daughter + suffix + extension
            with open(name, "w") as f:
                yaml.dump(d, f, default_flow_style=False)


    # Signal coupling yamls

    elif alp:

        center = -12.5
        width = np.linspace(0., 15., 10).tolist() + np.linspace(2., 6., 10).tolist()

        for i, w in enumerate(width):
            a = center - 0.5 * w
            b = center + 0.5 * w
            suffix = "_{}".format(i)
            d = alter_couplings(tempfile, 10.**a, 10.**b, suffix)

            name = daughter + '/' + daughter + suffix + extension
            with open(name, "w") as f:
                yaml.dump(d, f, default_flow_style=False)

    # Delete tempfile
    os.remove(tempfile)

    # Done!
