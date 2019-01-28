"""
Master module for all ColliderBit-related routines.
"""

import datetime

from files import *


def new_colliderbit_model(cb_output_dir, model):

    # Make sure that the output directories are in place
    dirs = [cb_output_dir + "/include/gambit/ColliderBit/models/", cb_output_dir + "/src/models/"]
    for d in dirs: mkdir_if_absent(d)

    for i, extension in enumerate([".hpp", ".cpp"]):
        # Take the template, replace the DATETIME and MODEL keys, and save the result
        old = os.getcwd() + "/Templates/ColliderBit_model" + extension
        new = dirs[i] + model + extension
        with open(old) as f_old, open(new, 'w') as f_new:
            for line in f_old:
                newline = re.sub("@MODEL@", model, line)
                newline = re.sub("@DATETIME@", datetime.datetime.now().strftime("%I:%M%p on %B %d, %Y"), newline)
                f_new.write(newline)
