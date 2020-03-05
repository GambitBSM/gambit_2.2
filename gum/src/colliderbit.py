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

    # Sort the higgses in numerical order - with the neutral ones first
    entry = sorted(entry, key=str.swapcase)
    listhiggses = ','.join(entry)

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

def get_higgs_invisibles(higgses, spheno_decays, particles, gambit_pdgs,
                         charged_higgses):
    """
    Generate a helper function "get_invisibles" for a new model. 
    Firstly, create a dictionary with potentially invisible decay products of 
    all Higgses, by PDG:

    { 25 : [ [12, -12], [14, -14], [16, -16], [darkmatter, darkmatter.conj] ] ,
      35 : [ [...] ] },

    then use it to create C++ code.

    TODO : only gets one LNP, does not know about multiple dark sectors
    TODO : assumes that the LNP is stable -- it might itself decay to visibles
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

    invis = {}
    all_parts = []

    for higgs in higgses:

        products = []

        # Get all decay products
        for decay in spheno_decays[higgs]:

            # Some quantities to query in a sec
            col = [int(colors[int(x)]) for x in decay]
            cha = [int(charges[int(x)]) for x in decay]
            stm = [int(sm[int(x)]) for x in decay]
            abspdgs = [abs(int(x)) for x in decay]

            # Save the products if they don't exist already
            for a in abspdgs:
                if not a in all_parts and sm[a] == False: 
                    all_parts.append(a)

            # Which of these are invisible? (And not BR to other Higgses)
            # Only want to save h -> LSP conj(LSP) too
            if (all([x == 1 for x in col]) and 
                all([x == 0 for x in cha]) and 
                all([x == 0 for x in stm]) and
                all([int(x) not in higgses for x in decay]) and
                abspdgs.count(abspdgs[0]) == len(decay)):

               products.append(decay)

        invis[higgs] = products

    # We now have a dictionary of all the possible states. Create a list of all
    # possible final states -- for *all* neutral Higgses, CP even and odd
    invisibles = list(set([tuple(sorted(x)) for y in invis.values() for x in y]))
    absinvs = [abs(int(i[0])) for i in invisibles]

    # "Visible" particles are those who are not invisible or Higgses, here...
    visibles = list( set(all_parts) - set(absinvs) - set(higgses) - 
                     set([abs(int(x)) for x in charged_higgses] ) )
    
    towrite = (
            "\n"
            "/// Helper function to work out if the LSP is invisible, "
            "and if so, which particle it is.\n"
            "std::vector<std::pair<str,str>> get_invisibles(const "
            "SubSpectrum& spec)\n"
            "{{\n"
            "// Get the lightest invisible particle (and antiparticle).\n"
            "std::pair<str,str> lnp = std::make_pair(\"{0}\", \"{1}\");\n"
            "double lnpmass = spec.get(Par::Pole_Mass, \"{2}\");\n"
    ).format(pdg_to_particle(abs(int(invisibles[0][0])), gambit_pdgs),
             pdg_to_particle(int(invisibles[0][0]), gambit_pdgs), 
             pdg_to_particle(int(invisibles[0][1]), gambit_pdgs))

    # Go through invisible particles -- find the lightest
    for pdg in invisibles[1:]:
        towrite += (
                "if (spec.get(Par::Pole_Mass, \"{0}\") < lnpmass)\n"
                "{{\n"
                "lnp = std::make_pair(\"{1}\", \"{2}\");\n"
                "lnpmass = spec.get(Par::Pole_Mass, \"{0}\");\n"
                "}}\n"
        ).format(pdg_to_particle(abs(int(pdg[0])), gambit_pdgs),
                 pdg_to_particle(int(pdg[0]), gambit_pdgs),
                 pdg_to_particle(int(pdg[1]), gambit_pdgs))

    towrite += (
            "// Work out if the lightest invisible particle is the LSP.\n"
            "bool inv_lsp = spec.get(Par::Pole_Mass, \"{0}\") > lnpmass "
    ).format(pdg_to_particle(visibles[0], gambit_pdgs))

    # TODO does not assume any mass ordering in eigenstates, not sure how
    # generic it is, so include them all atm.
    for pdg in visibles[1:]:
        towrite += (
                "and\n               spec.get(Par::Pole_Mass, \"{0}\") "
                "> lnpmass "
        ).format(pdg_to_particle(abs(pdg), gambit_pdgs))

    towrite += (
            ";\n"
            "// Check decays of at least one neutral higgs to it are "
            "kinematically possible.\n"
            "inv_lsp = (spec.get(Par::Pole_Mass, \"{0}\") > 2.*lnpmass"
    ).format(pdg_to_particle(higgses[0], gambit_pdgs))

    for higgs in higgses[1:]:
        towrite += (
                " or\n           spec.get(Par::Pole_Mass, \"{0}\") "
                "> 2.*lnpmass"
    ).format(pdg_to_particle(higgs, gambit_pdgs))

    towrite += (
            ");\n"
            "\n"
            "// Create a vector containing all invisible products "
            "of higgs decays.\n"
            "if (inv_lsp) return initVector<std::pair<str,str>>"
            "(lnp);\n"
            "// Return nothing if there's no contribution to h_inv.\n"
            "return std::vector<std::pair<str,str>>();\n"
            "}\n"
    )

    return indent(towrite)