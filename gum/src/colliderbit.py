"""
Master module for all ColliderBit-related routines.
"""

import datetime

from files import *
from setup import *


def new_colliderbit_model(cb_output_dir, model):
    """
    Creates a new ColliderBit model from the templates generated previously.
    """

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


def new_hct_switch(model_name, spectrum, neutral_higgses, gambit_pdgs):
    """
    Adds a new ModelInUse switch to the HiggsCouplingsTable
    routines in ColliderBit/src/ColliderBit_Higgs.cpp.
    """

    # Get the names of all neutral Higgses
    entry = []
    for higgs in neutral_higgses:
        entry.append("\""+pdg_to_particle(higgs, gambit_pdgs)+"\"")

    print neutral_higgses
    print entry


    # Sort the higgses in numerical order - with the neutral ones first
    entry = sorted(entry, key=str.swapcase)
    listhiggses = ','.join(entry)

    print entry
    print listhiggses

    towrite_src = (
                "else if (ModelInUse(\"{0}\"))\n"
                "{{\n"
                "  spectrum_dependency = &Dep::{1};\n"
                "  Higgses = initVector<str>({2});\n"
                "}}\n"
    ).format(model_name, spectrum, listhiggses)

    towrite_head = (
                 "    MODEL_CONDITIONAL_DEPENDENCY({0}, Spectrum, {1})\n"
    ).format(spectrum, model_name)

    return dumb_indent(6, towrite_src), towrite_head
