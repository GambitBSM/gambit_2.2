"""
Master module for all FeynRules related routines.
"""

from setup import *

def fr_part_to_gum_part(fr_bsm):
    """
    Converts all FeynRules::FRParticle to GUM::Particle.
    """
    
    bsm_list = []
    
    for i in xrange(len(fr_bsm)):
        part = fr_bsm[i]            
        bsm_list.append(Particle(part.spinX2(), part.pdg(), part.SC()))
    
    return bsm_list
