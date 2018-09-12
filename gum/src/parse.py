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

    def __init__(self, model_name, dm_candidate, mathpackage, lagrangian = None,
                 restriction = None, use_existing_spectrum = [False, None],
                 parent=None, children=None, friends=None):

        self.name = model_name
        self.dm_pdg = dm_candidate
        self.math = mathpackage
        self.restriction = None
        self.parent = None
        self.children = None
        self.friends = None
        self.LTot = lagrangian

        if use_existing_spectrum[0] == False:
            self.new_spectrum = True
            self.spec = "{0}_spectrum".format(model_name)
        else:
            self.new_spectrum = False
            self.spec = use_existing_spectrum[1]

        if restriction:
            self.restriction = restriction
        if parent:
            self.parent = parent
        if children:
            self.children = children
        if friends:
            self.friends = friends


    def is_orphan(self):
        """
        Check if model is an orphan i.e. it has no parent.
        """

        if parent:
            return True
        else:
            return False

    def has_friends(self):
        """
        Check if a model has friends or not.
        """

        if friends:
            return True
        else:
            return False

    def has_kids(self):
        """
        Check if a model has any children.
        """

        if children:
            return True
        else:
            return False

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

        if not 'dm_candidate' in data:
            raise GumError(("\n\nNo dark matter candidate specified. "
                            "Please check your .gum file."))

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
    dm_candidate = data['dm_candidate']

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

    # Model hierarchy stuff.
    restriction = None
    parent = children = friends = tf_p = tf_c = tf_f = None
    old_spectrum = spectrum_name = None
    if 'parent' in data:
        if not 'name' in data['parent']:
            raise GumError(("\n\nNo name given for parent function, please "
                            "check your .gum file."))
        parent = data['parent']['name']
        if not 'tf' in data['parent']:
            raise GumError(("\n\nNo translation function given for parent "
                            "function. Please check your .gum file."))
        tf_p = True

        # Assume we can use a parent's spectrum, unless explicitly specified
        old_spectrum = True
    if 'restriction' in math and mathpackage == 'feynrules':
        restriction = math['restriction']
    if 'children' in data:
        children = data['children']
    if 'friends' in data:
        friends = data['friends']
    if 'use_existing_spectrum' in data:
        old_spectrum = True
        spectrum_name = data['use_existing_spectrum']
        if not spectrum_name.endswith('_spectrum'):
            raise GumError(("\n\nExisting spectrum entry must end with "
                            "_spectrum. Please check your .gum file."))
    elif 'use_existing_spectrum' not in data:
        old_spectrum = False
        spectrum_name = "{0}_spectrum".format(gambit_model)

    dm_candidate = data['dm_candidate']

    gum_info = Inputs(gambit_model, dm_candidate, mathpackage, lagrangian,
                      restriction,
                      [old_spectrum, spectrum_name],
                      parent, children, friends)


    return gum_info, outputs

