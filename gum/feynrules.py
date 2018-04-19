"""
Master module for all FeynRules related routines.
"""

import os

from setup import *
from mathematica_interface import *

def check_fr_model(model_file):
    """
    Performs simple diagnostic routines to check the FeynRules file will
    function as expected.
    """

def load_fr_model(model_file):
    """
    Loads a FeynRules model.
    """

def obtain_feynrules_pdgs(model_file):
    """
    Obtains a dictionary of all BSM PDG codes.
    """

    feynrules_pdgs = {}
    return feynrules_pdgs

def obtain_feynrules_parameters(model_file):
    """
    Obtains a list of all new parameters in the FeynRules model.
    """
    
    feynrules_params = {}
    return feynrules_params

def parse_feynrules_model(model_file):
    """
    Parses a feynrules model file to extract all new parameters for use in
    GAMBIT. These are the parameters used to define the model in the Model
    Hierarchy.
    """

    pdgs = obtain_feynrules_pdgs(model_file)
    parameters = obtain_feynrules_parameters(model_file)

    return pdgs, parameters
  
def generate_chep(model_file):
    """
    Generates CalcHEP files from FeynRules.
    """
    
def generate_mg5(model_file):
    """
    Generate .ufo files from FeynRules.
    """

def generate_output(model_file, options=[]):
    """
    Runs the model through FeynRules, subject to certain options.
    """
    
    if 'chep' not in options:
        generate_chep(model_file)
    if 'mg5' not in options:
        generate_mg5(model_file)
    
