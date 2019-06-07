"""
Master module for all FeynRules related routines.
"""

from setup import *
import re 

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

def fr_params(paramlist, add_higgs):
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

    unsorted_params = []
    params = []

    # Add all parameters from the parameter list from FeynRules
    for i in xrange(len(paramlist)):
        p = paramlist[i]
        if (    (p.block() != 'YUKAWA')
            and (p.block() != 'SMINPUTS')
            and (p.block() != 'CKMBLOCK')):
            
            # Create a new instance of SpectrumParameter
            x = SpectrumParameter(p.name(), "dimensionless", block=p.block(), index=p.index())
            unsorted_params.append(x)

    # Now all of the parameters have been extracted, look to see if any of them
    # are elements of a matrix.

    # Firstly, get the name of each block
    blocks = list(set([i.block for i in unsorted_params]))

    # For each block, add the parameters to a dictionary...
    from collections import defaultdict
    blockdict = defaultdict(list)

    for i in blocks:
        for j in unsorted_params:
            if j.block == i:
                blockdict[i].append(j.name)

    # Parameters that GUM decides are (square) matrices
    matrices = {}

    # Go through each entry in the dictionary, and try and group them
    for k, v in blockdict.iteritems():

        # FeynRules splits a matrix M into M1x1, M1x2, ..., Mdxd.
        # Find any matches to M1x1 to start with.
        r = re.compile(r'(.*)1x1')
        first_entries = list(filter(r.match, v))

        # Now go through the list, to see how big the (square) matrix is.
        # Could generalise this to non-square matrices, if they're ever needed
        for i in first_entries:
            size = 1
            while (i[:-1] + str(size+1)) in v:
                size += 1
            matrices[i[:-3]] = "m{0}x{0}".format(str(size))

    # Delete duplicates from the original set of parameters
    keys = matrices.keys()
    added = [] # Dirty check to see a parameter's been added 
    for i in unsorted_params:
        present = False
        for k in keys:
            if k in i.name:
                present = True
                if k not in added:
                    added.append(k)
                    params.append(SpectrumParameter(k, "dimensionless", 
                                                    i.block, i.index, 
                                                    shape=matrices[k]))
            
        # If the parameter name doesn't match any of the matrix keys,
        # then just copy the parameter over from the unsorted list
        if not present:
            params.append(i)
        else: continue

    # Now add some Standard Model stuff that's in every SimpleSpectrum, for now.
    if add_higgs:
        params.append(SpectrumParameter("vev", "mass1", shape="scalar", sm=True))

    # Add gauge couplings and Yukawas here? TODO: check! 
    params.append(SpectrumParameter("g1", "dimensionless", shape="scalar", sm=True))
    params.append(SpectrumParameter("g2", "dimensionless", shape="scalar", sm=True))
    params.append(SpectrumParameter("g3", "dimensionless", shape="scalar", sm=True))
    params.append(SpectrumParameter("sinW2", "dimensionless", shape="scalar", sm=True))
    params.append(SpectrumParameter("Yd", "dimensionless", shape="m3x3", sm=True))
    params.append(SpectrumParameter("Yu", "dimensionless", shape="m3x3", sm=True))
    params.append(SpectrumParameter("Ye", "dimensionless", shape="m3x3", sm=True))
    
    return params