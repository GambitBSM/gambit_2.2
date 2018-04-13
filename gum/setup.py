"""
Master module containing class information, auto-writing codes, etc.
"""

import datetime
import os
from files import *

class GumError(Exception):
  pass

class Particle:
    """
    Particle class for internal use in GUM.
    """

    def __init__(self, chargex3, spinx2, pdg_code, own_conjugate):

        self.chargeX3 = chargex3
        self.spinX2 = spinx2
        self.PDG_code = pdg_code
        self.own_conjugate = own_conjugate

        # If a particle is self-conjugate -> the conjugate PDG code is same.
        if own_conjugate:
            self.Conjugate.PDG_code = self.PDG_code
            self.Conjugate.chargeX3 = self.chargeX3
        else:
            self.Conjugate.PDG_code = -self.PDG_code
            self.Conjugate.chargeX3 = -self.chargeX3

    def is_sc(self):
        """
        Boolean to flag if a particle is self conjugate.
        """

        return self.own_conjugate

    class Conjugate:
        """
        Conjugate sub-class for Particle class, e.g. DM.Conjugate.PDG_code
        for generalisation.
        """

        PDG_code = None
        chargex3 = None

        def __init__(self):
            self.PDG_code = self.PDG_code
            self.chargeX3 = self.chargeX3

def pdg_to_particle(pdg_code, pdg_dict):
    """
    Returns the particle name from the PDG code, from either
    a GAMBIT or CalcHEP dict, wrapped in quoatation marks.
    """
    

    for name, pdg_val in pdg_dict.iteritems():
        if pdg_code == pdg_val:
            return name
            
    # If not found -> throw error; gum doesn't know what to do.
    raise GumError(("No entry for PDG code " + str(pdg_code) + 
                    " in dictionary. Please check "
                    "gambit/config/particle_database.yaml "
                    "or your Mathematica file."))

class Vertex:
    """
    Vertex class for internal use in GUM.
    """

    def __init__(self):
        self.particles = []
        self.SM = True

    def is_sm(self):
        """
        Boolean to flag if an interaction is Standard Model.
        """
        return self.SM

    def num_particles(self):
        """
        Returns number of particles in vertex.
        """
        return len(self.particles)

class SpectrumParameter:
    """
    Spectrum Parameter class for use in GUM. Contains information about the
    Par::tag of parameter, read from a FeynRules or SARAH file, to be used
    in SpecBit and in Models.
    """

    def __init__(self, name, tag, shape=None, fullname=None,
                 sm=False, gb_input=None):
        self.name = name
        self.tag = tag
        self.shape = shape
        self.sm = sm
        if not fullname:
          self.fullname = name
        else:
          self.fullname = fullname
        if not gb_input:
          self.gb_in = name
        else:
          self.gb_in = gb_input
          
class BackendReq:
    """
    Backend requirement class for rollcall files. Assumes default
    behaviour is a BE_VARIABLE requiring no arguments. 
    """
    
    def __init__(self, capability, cpptype, tags=[], variable=True, 
                 args=[]):
      self.capability = capability
      self.cpptype = cpptype
      self.args = args
      self.tags = tags
      
      if variable:
        self.var = True
        self.args = None
      elif not variable:
        self.var = False
        self.args = args
        

class Dependency:
    """
    Dependency class used for rollcall files.
    """
    
    def __init__(self, name, cpptype):
      self.name = name
      self.cpptype = cpptype
      

def blame_gum(message):
    """
    Writes function to dump at the beginning of a new GAMBIT file.
    Blames GUM. Takes a message to describe the new file.
    """

    towrite = (
            "//   GAMBIT: Global and Modular BSM Inference Tool\n"
            "//   *********************************************\n"
            "///  \\file\n"
            "///\n"
            + message +
            "\n"
            "///\n"
            "///  Authors (add name and date if you modify):    \n"
            "///       *** Automatically created by GUM ***     \n"
            "///                                                \n"
            "///  \\author The GAMBIT Collaboration             \n"
            "///  \date " +
            datetime.datetime.now().strftime("%I:%M%p on %B %d, %Y") +
            "\n"
            "///                                                \n"
            "///  ********************************************* \n"
            "\n"
    )

    return towrite

def dumb_indent(numspaces, text):
   """
   Indents input text by specified number of spaces.
   """

   lines = text.splitlines(True)
   out = ""
   for line in lines:
     for i in range(0, numspaces):
       out += " "
     out += line

   return out

def indent(text, numspaces=0):
   """
   Increases indents of text every time it detects an open
   curly brace {, and decreases every time it detects a closed
   curly brace }. Does nothing if there is a {}.
   """

   lines = text.splitlines(True)
   out = ""

   for line in lines:


     if numspaces < 0:
       raise GumError(("Tried to indent a negative number of spaces."
                       "Please check for rogue braces."))

     num_open = line.count("{")
     num_closed = line.count("}")
     if num_open != num_closed:
       numspaces -= (num_closed*2)
     for i in range(0, numspaces):
        out += " "
     if num_open != num_closed:
       numspaces += (num_open*2)
     out += line

   return out
   
def reformat(location):
  """
  Reformat all text in a file (indents, etc.)
  """
  
