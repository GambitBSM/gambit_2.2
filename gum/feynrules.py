"""
Master module for all FeynRules related routines.
"""

## TO DO: WSTP/MathLink/etc. protocol. For now, just
## calling wolfram -script <model.m>

import os
from setup import *

def check_fr_model(model):
  """
  Performs simple diagnostic routines to check the FeynRules file will
  function as expected.
  """

def run_feynrules(script):
  """
  Runs FeynRules, via terminal, for a given .m file.
  """

  command = "wolfram -script " + script
  os.system(command)

def obtain_feynrules_pdgs(feynrules_path):
  """
  Obtains a dictionary of all BSM PDG codes.
  """

  feynrules_pdgs = {}
  return feynrules_pdgs

def obtain_feynrules_parameters(feynrules_path):
  """
  Obtains a list of all new parameters in the FeynRules model.
  """

  # Find where parameters begin to be defined in the model
  param_begin = open(feynrules_path, 'r').read().find('M$Parameters')

  if param_begin == -1:
    raise GumError("M$Parameters not found in FeynRules file. Please check.")

  with open(feynrules_path, 'r') as f:
    f.seek(param_begin)
    a = f.readlines()

  model_specific_parameters = []
  return model_specific_parameters

def parse_feynrules_model(feynrules_path):
  """
  Parses a feynrules model file to extract all new parameters for use in
  GAMBIT. These are the parameters used to define the model in the Model
  Hierarchy.
  """

  feynrules_pdgs = obtain_feynrules_pdgs(feynrules_path)
  model_specific_parameters = obtain_feynrules_parameters(feynrules_path)

  return feynrules_pdgs, model_specific_parameters
