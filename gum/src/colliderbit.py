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

def get_higgs_invisibles(higgses, spheno_decays, particles):
    """
    Returns a dictionary with potentially invisible decay products of 
    all Higgses, by PDG:

    { 25 : [ [12, -12], [14, -14], [16, -16], [darkmatter, darkmatter.conj] ] ,
      35 : [ [...] ] }

    which can be used by HiggsBounds routines.

    TODO : this should only be for neutral Higgses?
    """

    colors = {}
    charges = {}
    sm = {}
    for p in particles:
        colors[p.pdg()] = p.color()
        charges[p.pdg()] = p.chargeX3()
        charges[-p.pdg()] = -1*p.chargeX3()
        # Add antiparticles if distinct
        if p.color() != 1:
            colors[-p.pdg()] = -1*p.color()
        else:
            colors[-p.pdg()] = p.color()
        sm[p.pdg()] = p.SM()
        sm[-p.pdg()] = p.SM()

    d = {}

    for higgs in higgses:

        products = []

        # Get all decay products
        for decay in spheno_decays[higgs]:

            col = [int(colors[int(x)]) for x in decay]
            cha = [int(charges[int(x)]) for x in decay]
            stm = [int(sm[int(x)]) for x in decay]
            
            # Which of these are invisible? (And not BR to other Higgses)
            if (all([x == 1 for x in col]) and 
                all([x == 0 for x in cha]) and 
                all([x == 0 for x in stm]) and
                all([int(x) not in higgses for x in decay])):

               products.append(decay)

        d[higgs] = products

    return d