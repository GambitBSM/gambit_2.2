"""
Contains all routines for parsing input .gum file.
"""

import yaml

class Inputs:

  def __init__(self, model_name):
    self.name = model_name

def check_gum_file(inputfile):
  """
  Checks the input .GUM file for all necessary inputs.
  """
    
  print("Attempting to parse {0}...").format(inputfile)

  if inputfile.endswith(".gum"):
    pass

  else:
    raise GumError("Input filetype must be .gum.")

  with open(inputfile, "r") as f:
    try:
      data = yaml.load(f)
    except yaml.YAMLerror as exc:
      print(exc)
            
    if not 'new_model_name' in data:
    # This will just be a FR/SARAH *file* -> fine
      raise GumError(("No FeynRules script specified. "
                      "Please check your .gum file."))

    if not 'fr_script' in data:
    # This will just be a FR/SARAH *file* -> fine
      raise GumError(("No FeynRules script specified. "
                      "Please check your .gum file."))

    if not 'calchep_model' in data:
    # This will be linked to FR/SARAH directory -> fine
      raise GumError(("No CalcHEP model file directory specified. "
                      "Please check your .gum file."))

    if not 'dm_candidate' in data:
    # This will need to be specified in the .gum file.
      raise GumError(("No dark matter candidate specified. "
                      "Please check your .gum file."))

                            
  print("Parse successful.")

  return data
