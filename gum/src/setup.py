"""
Master module containing class information, auto-writing codes, etc.
"""

import datetime
import os

class GumError(Exception):
  pass

class Particle:
    """
    Particle class for internal use in GUM.
    """

    def __init__(self, name, antiname, spinx2, pdg_code, mass_name):

        self.name = name
        self.antiname = antiname
        self.spinX2 = spinx2
        self.PDG_code = pdg_code
        self.mass = mass_name

        self.own_conjugate = None
        if (name == antiname):
            self.own_conjugate = True
        else:
            self.own_conjugate = False

        # If a particle is self-conjugate -> the conjugate PDG code is same.
        if self.own_conjugate:
            self.Conjugate.PDG_code = self.PDG_code
        else:
            self.Conjugate.PDG_code = -self.PDG_code

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

        def __init__(self):
            self.PDG_code = self.PDG_code

def pdg_to_particle(pdg_code, pdg_dict):
    """
    Returns the particle name from the PDG code, from either
    a GAMBIT or CalcHEP dict, wrapped in quoatation marks.
    """
    
    for name, pdg_val in pdg_dict.iteritems():
        if pdg_code == pdg_val:
            return name
        
    """ 
    # If not found -> throw error; gum doesn't know what to do.
    raise GumError(("\n\nNo entry for PDG code " + str(pdg_code) + 
                    " in dictionary. Please check "
                    "gambit/config/particle_database.yaml "
                    "or your Mathematica file."))
    """
    # If not found -> return None & deal with this on case-by-case
    return None

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
      
