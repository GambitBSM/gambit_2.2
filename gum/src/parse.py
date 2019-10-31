"""
Contains all routines for parsing input .gum file.
"""

# TODO add mass dimension for new parameters in .gum file

import yaml
import re

from setup import *

class Inputs:
    """
    All the inputs from the .GUM file. Returns the master
    "gum" object used internally.
    """

    def __init__(self, model_name, base_model, mathpackage,
                 dm_candidate, mathname = None,
                 lagrangian = None, restriction = None):

        self.name = model_name.replace('-','_')
        self.base_model = base_model
        self.dm_pdg = dm_candidate
        self.math = mathpackage
        self.restriction = None
        self.LTot = lagrangian
        self.spec = "{0}_spectrum".format(model_name.replace('-','_'))

        # If we want the new GAMBIT model to have a different 
        # name than the model file from Mathematica
        if mathname:
            self.mathname = mathname
        else:
            self.mathname = model_name

        if restriction:
            self.restriction = restriction
        else:
            self.restriction = ''


class Outputs:
    """
    Outputs for GUM to write.
    """

    def __init__(self, mathpackage, calchep = False, pythia = False,
                 micromegas = False, spheno = False,
                 vevacious = False, ufo = False, collider_processes = None, 
                 multiparticles = None, pythia_groups = None):

        self.ch = calchep
        self.pythia = pythia
        self.mo = micromegas
        self.spheno = spheno
        self.vev = vevacious
        self.ufo = ufo
        self.collider_processes = collider_processes
        self.multiparticles = multiparticles
        self.pythia_groups = pythia_groups

        # Overwrite these, as the output does not exist.
        if mathpackage == 'feynrules':
            self.spheno = False
            self.vev = False

        # If Pythia is set then we also have UFO files, of course
        if pythia == True: self.ufo = True

    def bes(self):
        backends = []
        if self.ch: backends.append('calchep')
        if self.pythia: backends.append('pythia')
        if self.mo: backends.append('micromegas')
        if self.spheno: backends.append('spheno')
        if self.vev: backends.append('vevacious')
        if self.ufo: backends.append('ufo')
        return backends


def check_gum_file(inputfile):
    """
    Checks the input .GUM file for all necessary inputs.
    """

    print("Attempting to parse {0}...").format(inputfile)

    if inputfile.endswith(".gum"):
        pass
    else:
        raise GumError("\n\nInput filetype must be .gum.")

    with open(inputfile, "r") as f:
        try:
            data = yaml.load(f)
        except yaml.YAMLerror as exc:
            print(exc)

        if not 'math' in data:
            raise GumError(("\n\n'math' node needed in .gum file."))

        if not 'package' in data['math']:
        # Don't know what to run...!
            raise GumError(("\n\nNo mathpackage input - what do you expect "
                            "GUM to do? Please check your .gum file. "
                            "Supported entries: sarah, feynrules."))

        if data['math']['package'] not in ["sarah", "feynrules"]:
            raise GumError(("\n\nYou must specify which mathpackage you want "
                            "GUM to use. Please check your .gum file. "
                            "Supported entries: sarah, feynrules."))

        if not 'model' in data['math']:
            raise GumError(("\n\nNo model file specified. "
                            "Please check your .gum file."))


    print("Parse successful.")

    return data

def fill_gum_object(data):
    """
    Returns a model of type Inputs for GUM to work with. 'data' is the
    parsed data from check_gum_file.
    """

    math = data['math']
    if 'lagrangian' in data['math']:
        # Check the Lagrangian makes sense (i.e. is all alphanumeric)
        lagrangian = data['math']['lagrangian']
        L = lagrangian.split('+')
        for l in L:
            if not l.strip(' ').isalnum():
                raise GumError(("Non-alphanumeric character detected in "
                                " the Lagrangian. Please check your .gum "
                                "file."))
    else:
        lagrangian = "LTotal"
    mathpackage = math['package']

    gambit_model = math['model']

    # Overwrite the GAMBIT model if specified
    mathname = ""
    if 'gambit_opts' in data:
        if 'model_name' in data['gambit_opts']:
            mathname = gambit_model
            gambit_model = data['gambit_opts']['model_name']

    # FeynRules specific -- a "base" model to build a pheno model on top of.
    # Tyically this is the SM, plus the BSM contribution defined in a separate file.
    if 'base_model' in data['math']:
        base_model = data['math']['base_model']
    else:
        base_model = ""

    if 'dm_candidate' in data:
        dm_candidate = data['dm_candidate']
    else:
        dm_candidate = None


    backends = ['calchep', 'pythia', 'spheno', 'ufo',
                'micromegas', 'vevacious']

    opts = {}
    # The outputs GUM should hook up to GAMBIT, if specified
    if 'output' in data:
        for i in backends:
            if i in data['output']:
                opts[i] = data['output'][i]
            else:
                opts[i] = False
        if 'collider_processes' in data['output']:
            opts['collider_processes'] = data['output']['collider_processes']
        if 'multiparticles' in data['output']:
            opts['multiparticles'] = data['output']['multiparticles']
        if 'pythia_groups' in data['output']:
            opts['pythia_groups'] = data['output']['pythia_groups']

    # Default: if unspecified, write everything.
    else:
        for i in backends:
            opts[i] = True
    
    outputs = Outputs(mathpackage, **opts)

    # If the user wants MicrOMEGAs output but hasn't specified a DM candidate
    if not dm_candidate and outputs.mo:
        raise GumError(("\n\nYou have asked for MicrOMEGAs output but have not "
                        "specified which particle is meant to be the DM "
                        "candidate! Please add an entry to your .gum file "
                        "like:\n\ndm_candidate: 9900001 # <--- insert the "
                        "desired PDG code here!!\n")) 


    # FeynRules restriction files
    restriction = None
    if 'restriction' in math and mathpackage == 'feynrules':
        restriction = math['restriction']

    gum_info = Inputs(gambit_model, base_model, mathpackage, 
                      dm_candidate, mathname, lagrangian, restriction)


    return gum_info, outputs

