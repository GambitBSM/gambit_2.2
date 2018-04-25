"""
Master module for all SARAH related routines.
"""

import os

from setup import *
from mathematica_interface import *

def load_sarah_model(model_file):
    """
    Loads a SARAH model.
    """
    
def check_sarah_model(model_file):
    """
    Performs simple diagnostic routines to check the SARAH file will
    function as expected.
    """

def obtain_sarah_pdgs(model_file):
    """
    Obtains a dictionary of all BSM PDG codes.
    """

    sarah_pdgs = {}
    return sarah_pdgs

def obtain_sarah_parameters(model_file):
    """
    Obtains a list of all new parameters in the SARAH model.
    """
    
    sarah_params = {}
    return sarah_params

def parse_sarah_model(model_file):
    """
    Parses a sarah model file to extract all new parameters for use in
    GAMBIT. These are the parameters used to define the model in the Model
    Hierarchy.
    """

    pdgs = obtain_sarah_pdgs(model_file)
    parameters = obtain_sarah_parameters(model_file)

    return pdgs, parameters
  
def generate_chep(model_file):
    """
    Generates CalcHEP files from SARAH.
    """
    
def generate_mg5(model_file):
    """
    Generate .ufo files from SARAH.
    """

def generate_spheno(model_file):
    """
    Generate SPheno output from SARAH.
    """
    
def generate_vevacious(model_file):
    """
    Generate VEVacious output from SARAH.
    """
    
def generate_fs(model_file):
    """
    Generate FlexibleSUSY output from SARAH.
    """

def generate_output(model_file, options=[]):
    """
    Runs the model through SARAH, subject to certain options.
    """
    
    if 'chep' not in options:
        generate_chep(model_file)
    if 'mg5' not in options:
        generate_mg5(model_file)
    if 'spheno' not in options:
        generate_spheno(model_file)
    if 'fs' not in options:
        generate_fs(model_file)    
    if 'vevacious' not in options:
        generate_vevacious(model_file)
