"""
Contains all routines for parsing input .gum file.
"""

import yaml
from setup import *

class Inputs:

    def __init__(self, model_name, dm_candidate, mathpackage, 
                 use_existing_spectrum = [False, None], parent=None, 
                 children=None, friends=None):
                 
        self.name = model_name
        self.dm_name = dm_candidate
        
        self.parent = None
        self.children = None
        self.friends = None
        
        if use_existing_spectrum[0] == False:
            self.new_spectrum = True
            self.spec = "{0}_spectrum".format(model_name)
        else:
            self.new_spectrum = False
            self.spec = use_existing_spectrum[1]
            
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
                             
        if not 'mathpackage' in data:
        # Don't know what to run...!
            raise GumError(("No mathpackage input - what do you expect "
                            "GUM to do? Please check your .gum file. "
                            "Supported entries: sarah, feynrules."))
                                    
        if data['mathpackage'] not in ["sarah", "feynrules"]:
            raise GumError(("You must specify which mathpackage you want "
                            "GUM to use. Please check your .gum file. "
                            "Supported entries: sarah, feynrules."))
  
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
  
def fill_gum_object(data):
    """
    Returns a model of type Inputs for GUM to work with. 'data' is the 
    parsed data from check_gum_file.
    """
    
    gambit_model = data['new_model_name']
    mathpackage = data['mathpackage']
    # Model hierarchy stuff.
    parent = children = friends = tf_p = tf_c = tf_f = None 
    old_spectrum = spectrum_name = None
    if 'parent' in data:
        if not 'name' in data['parent']:
            raise GumError(("No name given for parent function, please check "
                            "your .gum file."))
        parent = data['parent']['name']
        if not 'tf' in data['parent']:
            raise GumError(("No translation function given for parent function."
                            " Please check your .gum file."))
        tf_p = True
        
        # Assume we can use a parent's spectrum, unless explicitly specified
        old_spectrum = True
    if 'children' in data:
        children = data['children']
    if 'friends' in data:
        friends = data['friends']
    if 'use_existing_spectrum' in data:
        old_spectrum = True
        spectrum_name = data['use_existing_spectrum']
        if not spectrum_name.endswith('_spectrum'):
            raise GumError(("Existing spectrum entry must end with _spectrum. "
                            "Please check your .gum file."))
    elif 'use_existing_spectrum' not in data:
        old_spectrum = False
        spectrum_name = "{0}_spectrum".format(gambit_model)
      
    dm_candidate = data['dm_candidate']
    
    gum_info = Inputs(gambit_model, dm_candidate, mathpackage, [old_spectrum, spectrum_name],
                      parent, children, friends)

    
    return gum_info
    
