"""
Contains all routines for parsing input .gum file.
"""

import yaml

from setup import *

class Inputs:
    """
    All the inputs from the .GUM file. Returns the master
    "gum" object used internally.
    """

    def __init__(self, model_name, mathpackage,
                 dm_candidate,
                 lagrangian = None, restriction = None):

        self.name = model_name
        self.dm_pdg = dm_candidate
        self.math = mathpackage
        self.restriction = None
        self.LTot = lagrangian
        self.spec = "{0}_spectrum".format(model_name)

        if restriction:
            self.restriction = restriction
        else:
            self.restriction = ''


class Outputs:
    """
    Outputs for GUM to write.
    """

    def __init__(self, mathpackage, calchep = False, pythia = False,
                 micromegas = False, spheno = False, flexiblesusy = False,
                 vevacious = False):

        self.ch = calchep
        self.ufo = pythia
        self.mo = micromegas
        self.sph = spheno
        self.fs = flexiblesusy
        self.vev = vevacious

        # Overwrite these, as the output does not exist.
        if mathpackage == 'feynrules':
            self.mo = False
            self.sph = False
            self.fs = False
            self.vev = False

    def bes(self):
        backends = []
        if self.ch: backends.append('calchep')
        if self.ufo: backends.append('pythia')
        if self.mo: backends.append('micromegas')
        if self.sph: backends.append('spheno')
        if self.fs: backends.append('flexiblesusy')
        if self.vev: backends.append('vevacious')
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
        lagrangian = data['math']['lagrangian']
    else:
        lagrangian = "LTotal"
    mathpackage = math['package']
    gambit_model = math['model']
    if 'dm_candidate' in data:
        dm_candidate = data['dm_candidate']
    else:
        dm_candidate = None

    backends = ['calchep', 'pythia', 'spheno', 'flexiblesusy',
                'micromegas', 'vevacious']

    opts = {}
    # The outputs GUM should hook up to GAMBIT, if specified
    if 'output' in data:
        for i in backends:
            if i in data['output']:
                opts[i] = data['output'][i]
            else:
                opts[i] = False

    # Default: if unspecified, write everything.
    else:
        for i in backends:
            opts[i] = True

    outputs = Outputs(mathpackage, **opts)

    # FeynRules restriction files
    restriction = None
    if 'restriction' in math and mathpackage == 'feynrules':
        restriction = math['restriction']

    gum_info = Inputs(gambit_model, mathpackage,  dm_candidate, lagrangian,
                      restriction)


    return gum_info, outputs

