"""
Master module for all SARAH related routines.
"""

from setup import *

def sarah_part_to_gum_part(sarah_bsm):
    """
    Converts all SARAH::SARAHParticle to GUM::Particle.
    """
    
    bsm_list = []
    
    for i in xrange(len(sarah_bsm)):
        part = sarah_bsm[i]            
        bsm_list.append(Particle(part.spinX2(), part.pdg(), 
                                 part.SC(), part.mass()))
    
    return bsm_list
