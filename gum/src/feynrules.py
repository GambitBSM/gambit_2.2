"""
Master module for all FeynRules related routines.
"""

from setup import *

def fr_part_to_gum_part(fr_bsm):
    """
    Converts all FeynRules::FRParticle to GUM::Particle.
    Returns whether or not to add a dependency on StandardModel_Higgs.
    """
    
    bsm_list = []
    add_higgs = False

    all_pdgs = [x.pdg() for x in fr_bsm]
    # If the SM higgs included but not a second Higgs, 
    # then we want a dependency on StandardModel_Higgs.
    # If it is a 2+HDM then we don't.
    if 25 in all_pdgs and 35 not in all_pdgs: 
        add_higgs = True

    for i in xrange(len(fr_bsm)):
        part = fr_bsm[i]
        bsm_list.append(Particle(part.name(), part.antiname(),
        						 part.spinX2(), part.pdg(), 
                                 part.mass()))
    
    return bsm_list, add_higgs

def fr_params(paramlist):
    """
    Removes all Standard Model parameters from those we wish
    to add to the GAMBIT model. This utilises the 'BlockName'
    tag (i.e. the SLHA Block the parameters will be placed in)
    to distinguish between Standard Model parameters, and 
    parameters the user has implemented in their own FeynRules
    model file.

    Behaviour mimics the SM.fr file that is delivered with 
    FeynRules by default. The following blocks will not be 
    considered BSM: 'SMINPUTS', 'YUKAWA' or 'CKMBLOCK'.
    Everything else is fair game. 

    Default assumption is that everything is a dimensionless
    parameter. 
    TODO. Since there are no spectra, this should be ok.
    Maybe worth adding a new parameter tag, like Par::unspecified.
    """

    params = []

    # Add all parameters from the parameter list from FeynRules
    for i in xrange(len(paramlist)):
        p = paramlist[i]
        if (    (p.block() != 'YUKAWA')
            and (p.block() != 'SMINPUTS')
            and (p.block() != 'CKMBLOCK')):
            
            # Create a new instance of SpectrumParameter
            x = SpectrumParameter(p.name(), "dimensionless")
            params.append(x)

    # Add gauge couplings and Yukawas here? Unsure. 
    # They're vetoed by the SMINPUTS filter above.
    """
    params.append(SpectrumParameter("g1", "dimensionless", "scalar", sm=True))
    params.append(SpectrumParameter("g2", "dimensionless", "scalar", sm=True))
    params.append(SpectrumParameter("g3", "dimensionless", "scalar", sm=True))
    params.append(SpectrumParameter("sinW2", "dimensionless", "scalar", sm=True))
    params.append(SpectrumParameter("Yd", "dimensionless", "m3x3", sm=True))
    params.append(SpectrumParameter("Yu", "dimensionless", "m3x3", sm=True))
    params.append(SpectrumParameter("Ye", "dimensionless", "m3x3", sm=True))
    """


    return params