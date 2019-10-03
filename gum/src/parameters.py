#!/usr/bin/env python
#
#  GUM: GAMBIT Universal Models
#  ****************************
#  \file
#
#  Master module for all routines relating 
#  to parameters, manipulation, etc.
#
#  *************************************
#
#  \author Sanjay Bloor
#          (sanjay.bloor12@imperial.ac.uk)
#  \date 2017, 2018, 2019 
#
#  \author Tomas Gonzalo
#          (tomas.gonzalo@monash.edu)
#  \date 2019 Aug
#
#  **************************************

from setup import *
import re

###############
## FEYNRULES ##
###############

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
            x = SpectrumParameter(p.name(), "dimensionless", block=p.block(), 
                                  index=p.index())
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
    params.append(SpectrumParameter("g1", "dimensionless", block="GAUGE", index=1, shape="scalar", sm=True))
    params.append(SpectrumParameter("g2", "dimensionless", block="GAUGE", index=2, shape="scalar", sm=True))
    params.append(SpectrumParameter("g3", "dimensionless", block="GAUGE", index=3, shape="scalar", sm=True))
    params.append(SpectrumParameter("sinW2", "dimensionless", shape="scalar", sm=True))
    params.append(SpectrumParameter("Yd", "dimensionless", block="YD", shape="m3x3", sm=True))
    params.append(SpectrumParameter("Yu", "dimensionless", block="YU", shape="m3x3", sm=True))
    params.append(SpectrumParameter("Ye", "dimensionless", block="YE", shape="m3x3", sm=True))
    
    return params

###########
## SARAH ##
###########

def sarah_part_to_gum_part(sarah_bsm):
    """
    Converts all SARAH::SARAHParticle to GUM::Particle.
    """
    
    bsm_list = []
    add_higgs = False

    all_pdgs = [x.pdg() for x in sarah_bsm]
    # If the SM higgs included but not a second Higgs, 
    # then we want a dependency on StandardModel_Higgs.
    # If it is a 2+HDM then we don't.
    if 25 in all_pdgs and 35 not in all_pdgs: 
        add_higgs = True
    
    for i in xrange(len(sarah_bsm)):
        part = sarah_bsm[i]            
        bsm_list.append(Particle(part.name(), part.antiname(),
        			 part.spinX2(), part.pdg(), 
                                 part.mass(), part.alt_name(),
                                 part.alt_mass()))
    
    return bsm_list, add_higgs

def sarah_params(paramlist, add_higgs):
    """
    Removes all Standard Model parameters from those we wish
    to add to the GAMBIT model. This utilises the 'BlockName'
    tag (i.e. the SLHA Block the parameters will be placed in)
    to distinguish between Standard Model parameters, and 
    parameters the user has implemented in their own FeynRules
    model file.

    Behaviour mimics the default imported SARAH conventions, which
    for some strange reason uses the description tag of parameters
    as protected variables. 
    The procudure is:
    - read in a new model
    - read in $SARAH/Models/parameters.m (and particles.m)
    - overwrite any parameter definitions that have the same
     *descriptions* in both files with that from the defaults.
    Anything that is given in the following blocks will be ignored:
    'SMINPUTS', 'CKMBLOCK', 'SM'
    Everything else is fair game. 

    Default assumption is that everything is a dimensionless
    parameter. 
    TODO. This is a problem for SARAH, and actual spectrum objects.
    Hmmm......
    Maybe worth adding a new parameter tag, like Par::unspecified.
    """

    unsorted_params = []
    params = []

    # Add all parameters from the parameter list from SARAH
    for i in xrange(len(paramlist)):
        p = paramlist[i]
        if (    (p.block() != 'SM')
            and (p.block() != 'SMINPUTS')
            and (p.block() != 'VCKM')
            #and (p.is_output() != False)) TODO: TG: This line was killing most of the parameters
            ):
            
            # Create a new instance of SpectrumParameter
            # TODO: dimensionless atm! 
            x = SpectrumParameter(p.name(), "dimensionless", block=p.block(),
                                  index=p.index(), alt_name = p.alt_name(),
                                  bcs = p.bcs(), shape = p.shape(), 
                                  is_output = p.is_output(), is_real = p.is_real())
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

    for par in unsorted_params:
      params.append(par)
    
    # Now add some Standard Model stuff that's in every SimpleSpectrum, for now.
    if add_higgs:
        params.append(SpectrumParameter("vev", "mass1", shape="scalar", sm=True, is_real=True))

    # Add gauge couplings and Yukawas here? TODO: check! 
    params.append(SpectrumParameter("g1", "dimensionless", block="GAUGE", index=1, shape="scalar", sm=True, is_real=True))
    params.append(SpectrumParameter("g2", "dimensionless", block="GAUGE", index=2, shape="scalar", sm=True, is_real=True))
    params.append(SpectrumParameter("g3", "dimensionless", block="GAUGE", index=3, shape="scalar", sm=True, is_real=True))
    params.append(SpectrumParameter("sinW2", "dimensionless", shape="scalar", sm=True, is_real=True))
    params.append(SpectrumParameter("Yd", "dimensionless", block="YD", shape="m3x3", sm=True, is_real=True))
    params.append(SpectrumParameter("Yu", "dimensionless", block="YU", shape="m3x3", sm=True, is_real=True))
    params.append(SpectrumParameter("Ye", "dimensionless", block="YE", shape="m3x3", sm=True, is_real=True))
    
    return params


##################
# Other routines #
##################

def sort_params_by_block(parameters):
    """
    Returns a dictionary of the LH blocks of a new model,
    with all entries.

    The dict looks like: 

    {
        # A block with many entries
        SMINPUTS : { 1: alphainv, 2: GF, 3: alphaS, ... }, 

        # A block with just one entry, i.e. matrices:
        YE : { matrix: 3x3 },         # e.g. Yukawas

        # A matrix block with a different outputname in SPheno,
        # e.g. for the THDM
        SCALARMIX : { mixingmatrix: 2x2, outputname: ZH }
        ...
    }
    """

    params_by_block = {}

    for par in parameters:

        # If there's no block, we can't really add it to this dict.
        if not par.block:
            continue

        # If the parameter 

        shape = par.shape

        matrix = False
        if shape:
            if shape.startswith('m'): 
                matrix = True

        # If it's a matrix then it will be a new block
        if matrix:
            if par.is_output:
                newentry = { "mixingmatrix": par.shape[1:], "outputname": par.name }
            else:
                newentry = { "matrix": par.shape[1:], "outputname": par.name }
            params_by_block[par.block] = newentry

        # If it's not a matrix and is a new block, then create the entry
        elif par.block not in params_by_block:
            newentry = { par.index : par.name }
            params_by_block[par.block] = newentry

        # If the block already exists, add to it
        elif par.block in params_by_block:
            params_by_block[par.block][par.index] = par.name

    return params_by_block

def parameter_reality_dict(parameters):
    """
    Get a dict of whether a parameter is most definitely real, or not.
    """

    reality_dict = {}

    for param in parameters:
        reality_dict[param.name] = param.is_real

    return reality_dict

def spheno_dependencies(sphenodeps):
    """
    Clean the spheno dependencies scraped from SARAH.
    Here 'clean' means tidying up the CForm output Mathematica gives us.
    Returns a dict of parameter-definition key-value pairs.
    """    

    deps = {}

    # each p is an instance of SARAHParameter
    for p in sphenodeps:

        name = p.name() # As known to SPheno
        # Don't need the block. index etc, we're using internal SPheno params here

        description = p.alt_name()

        # Common replacements to actual C(++) syntax
        # Probably not complete, but the main culprits *should* be here
        description = re.sub('ArcSin', 'asin', description)
        description = re.sub('ArcCos', 'acos', description)
        description = re.sub('ArcTan', 'atan', description)
        description = re.sub('Abs', 'abs', description) 
        description = re.sub('Sin', 'sin', description)
        description = re.sub('Cos', 'cos', description)
        description = re.sub('Tan', 'tan', description)
        #description = re.sub('Power', 'pow', description)
        description = re.sub('Sqrt', 'sqrt', description)

        # If there's something that looks like param(i,j) make this (param)(i,j)
        description = re.sub(r'([a-zA-Z]+)\(([0-9]),([0-9])\)',r'(*\1)(\2,\3)', description)

        # If there's Power(param, num) -> pow(*param, num)
        description = re.sub(r'Power\(([a-zA-Z]+),([0-9])\)',r'pow(*\1,\2)', description)

        deps[name] = description

    return deps