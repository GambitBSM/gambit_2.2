#!/usr/bin/env python
#
#  GUM: GAMBIT Universal Models
#  ****************************
#  \file
#
#  Contains all routines for parsing input .gum files and 
#  the SARAH/FeynRules model files.
#
#  *************************************
#
#  \author Sanjay Bloor
#          (sanjay.bloor12@imperial.ac.uk)
#  \date 2017, 2018, 2019, 2020
#
#  \author Pat Scott
#          (pat.scott@uq.edu.au)
#  \date 2018, 2019
#
#  \author Tomas Gonzalo
#          (tomas.gonzalo@monash.edu)
#  \date 2020 Jan
#
#  **************************************


# TODO add mass dimension for new parameters in .gum file

import yaml
import re

from setup import *

"""
.GUM FILE PARSING
"""

class Inputs:
    """
    All the inputs from the .GUM file. Returns the master
    "gum" object used internally.
    """

    def __init__(self, model_name, base_model, mathpackage,
                 wimp_candidate, mathname = None,
                 lagrangian = None, restriction = None):

        self.name = model_name.replace('-','_')
        self.base_model = base_model
        self.dm_pdg = wimp_candidate
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
                 vevacious = False, ufo = False, options = {}):
#                  , collider_processes = None, 
#                 multiparticles = None, pythia_groups = None):

        self.ch = calchep
        self.pythia = pythia
        self.mo = micromegas
        self.spheno = spheno
        self.vev = vevacious
        self.ufo = ufo
#        self.collider_processes = collider_processes
#        self.multiparticles = multiparticles
#        self.pythia_groups = pythia_groups
        self.options = options

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
        if inputfile.endswith(".mug"):
            raise GumError(("\n\nGUM called with a .mug file in normal mode --"
                            " you probably want to call gum with the -r flag:"
                            "\n\n  ./gum -r " + inputfile + "\n"))
        else:
            raise GumError("\n\nInput filetype must be .gum.")

    with open(inputfile, "r") as f:
        try:
            data = yaml.load(f)
        except yaml.YAMLError as exc:
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

        if not 'output' in data:
        # Don't know what to generate!
            raise GumError(("\n\nNo output specified! You need to tell GUM "
                            "what it is you'd like it to do!\n"
                            "Please change your .gum file!"))

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
    # Tyically this is the SM, plus the BSM contribution defined in a 
    # separate file.
    if 'base_model' in data['math']:
        base_model = data['math']['base_model']
    else:
        base_model = ""

    if 'wimp_candidate' in data:
        wimp_candidate = data['wimp_candidate']
    else:
        wimp_candidate = None

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
        if all(value == False for value in opts.values()):
            raise GumError(("\n\nAll backend output set to false in your .gum "
                            "file.\nGive GUM something to do!\n"
                            "Please change your .gum file."))

    options = {}
    # Options for the outputs declared
    if 'output_options' in data:
        for output in data['output_options']:
            if output not in opts.keys():
                raise GumError(("\n\nOptions given to output " + output + " "
                                "which is not declared as gum output.\n"
                                "Please change your .gum file."))
            options[output] = data['output_options'][output]

        #if 'collider_processes' in data['output']:
        #    opts['collider_processes'] = data['output']['collider_processes']
        #if 'multiparticles' in data['output']:
        #    opts['multiparticles'] = data['output']['multiparticles']
        #if 'pythia_groups' in data['output']:
        #    opts['pythia_groups'] = data['output']['pythia_groups']

    # If we've got this far, we'll also force some decays to be written,
    # either by SPheno or by CalcHEP.
    # N.B. vevacious is conditional on SPheno
    # This now means if any of: Pythia(ufo) or MicrOMEGAs are requested, 
    # then we'll default to activating CalcHEP (unless SPheno requested)
    set_calchep = True
    if mathpackage == 'sarah' and opts['spheno'] == True:
        set_calchep = False
    if set_calchep: opts['calchep'] = True
    
    outputs = Outputs(mathpackage, options=options, **opts)

    # If the user wants MicrOMEGAs output but hasn't specified a DM candidate
    if not wimp_candidate and outputs.mo:
        raise GumError(("\n\nYou have asked for MicrOMEGAs output but have not "
                        "specified which particle is meant to be the DM "
                        "candidate! Please add an entry to your .gum file "
                        "like:\n\nwimp_candidate: 9900001 # <--- insert the "
                        "desired PDG code here!!\n")) 

    # FeynRules restriction files
    restriction = None
    if 'restriction' in math and mathpackage == 'feynrules':
        restriction = math['restriction']

    gum_info = Inputs(gambit_model, base_model, mathpackage, 
                      wimp_candidate, mathname, lagrangian, restriction)


    return gum_info, outputs

"""
FEYNRULES PARSING
"""

def parse_feynrules_model_file(fr_file_path, outputs):
    """
    Parses a FeynRules model file. Checks for the following:
        - Every parameter has an LH block and an
        - No particles have an external name with underscores etc.
        - Every parameter has an interaction order specified *if* the user
          requests UFO output
        - ComplexParameters and CalcHEP output
    """

    payattn = False

    # Read the input in
    with open(fr_file_path, 'r') as f:
        lines = f.readlines()

    # Flatten the string
    contents = "".join(lines).replace("\n","")

    # Parse the contents

    # 1. Parameters
    blocks = {}
    interactionorders = {}

    numbraces = 0

    # Remove all the stuff before parameters
    s1 = contents[contents.find('M$Parameters'):]

    # Search through the string and count the number of curly braces. 
    # When we get to zero then we're done.
    numbraces = 0
    started = False
    commenting = False
    s2 = ""
    for char in s1:
        if numbraces > 0: started = True
        if char == "{" and not commenting:   numbraces += 1
        elif char == "}" and not commenting: numbraces -= 1
        # Don't count braces if we're in a comment, just in case someone is 
        # truly twisted
        elif char == "\"": commenting = not commenting
        if started: s2 += char
        if numbraces == 0 and started: break

    # Get the param names
    params = []
    for i in s2.split('==')[:-1]:
        params.append( i.strip(' ').split(' ')[-1])
    
    # Get the contents between each parameter
    matches = []
    for i in range(len(params)-1):
        pat = r'{}(.*?){}'.format(params[i], params[i+1])
        matches.append( re.search(pat, s2).group(1) )
    matches.append( s2[s2.find(params[-1]):])

    for i in range(len(matches)):
        match = matches[i]
        blockmatch = re.search(r'BlockName\s*->(.*?),', match)
        ordermatch = re.search(r'OrderBlock\s*->(.*?),', match)
        if blockmatch:
            if ordermatch: 
                blocks[ params[i] ] = { blockmatch.group(1) : 
                                        ordermatch.group(1) }
            else:
                blocks[ params[i] ] = blockmatch.group(1)
        # else:
        #     raise GumError(("No BlockName specified for the parameter "
        #                     "{0}. GUM and GAMBIT need this information "
        #                     "so please change your .fr file!"
        #                     ).format(params[i]))
        interactionmatch = re.search(r'InteractionOrder\s*->\s*\{(.*?),(.*?)\}',
                                     match)
        if interactionmatch:
            interactionorders[ params[i] ] = interactionmatch.group(1)
        else:
            interactionorders[ params[i] ] = "NULL"

    # Check for InteractionOrder
    # TODO need to specify *which* parameters actually need these...
    if outputs.ufo:
        for param, orders in interactionorders.iteritems():
            if orders == "NULL":
                raise GumError(("No interaction order specified for the "
                                "parameter {0}, and you have asked for UFO "
                                "output. If you want to use UFO output, every "
                                "parameter needs an InteractionOrder specified."
                                " Please consult the FeynRules and/or MadGraph "
                                "manuals for details.").format(param))

    # 2. Particles
    numbraces = 0

    # Remove all the stuff before parameters
    s3 = contents[contents.find('M$ClassesDescription'):]

    # Search through the string and count the number of curly braces. 
    # When we get to zero then we're done.
    numbraces = 0
    started = False
    commenting = False
    s4 = ""
    for char in s3:
        if numbraces > 0: started = True
        if char == "{" and not commenting:   numbraces += 1
        elif char == "}" and not commenting: numbraces -= 1
        # Don't count braces if we're in a comment, just in case someone is 
        # truly twisted
        elif char == "\"": commenting = not commenting
        if started: s4 += char
        if numbraces == 0 and started: break

    # Get the particles names
    particles = []
    for i in s4.split('==')[:-1]:
        particles.append( i.strip(' ').split(' ')[-1])
    
    # Get the contents between each parameter
    partmatches = []
    for i in range(len(particles)-1):
        p1 = particles[i].replace('[',r'\[').replace(']',r'\]')
        p2 = particles[i+1].replace('[',r'\[').replace(']',r'\]')
        pat = r'{}(.*?){}'.format(p1, p2)
        partmatches.append( re.search(pat, s4).group(1) )
    partmatches.append( s4[s4.find(particles[-1]):])

    for i in range(len(partmatches)):
        # PDG codes for ghosts don't matter
        if re.match(r'U\[\d*\]', particles[i]):
            continue
        match = partmatches[i]
        # Unphysical fieldsd don't matter either
        unphysmatch = re.search(r'Unphysical\s*->\s*True', match)
        if unphysmatch:
            continue
        partmatch = re.search(r'ClassName\s*->(.*?),', match)
        if not partmatch:
            raise GumError(("The particle with description {0} "
                            "does not have a class name. Your FeynRules "
                            "file shouldn't work..."
                             ).format(particles[i]))
        # Try to match a set of PDG codes first
        pdgmatch = re.search(r'PDG\s*->\s*\{(.*?)\}\s*', match)
        if not pdgmatch:
            pdgmatch = re.search(r'PDG\s*->(.*?)\s*,', match)

        if not pdgmatch:
            raise GumError(("Particle {0} does not have a PDG code; "
                            "please give it one in your .fr file, as "
                            "GUM and GAMBIT need this information."
                            ).format(particles[i]))
