#!/usr/bin/env python
#
#  GUM: GAMBIT Universal Models
#  ****************************
#  \file
#
#  Master module for all Models related routines.
#
#  *************************************
#
#  \author Sanjay Bloor
#          (sanjay.bloor12@imperial.ac.uk)
#  \date 2019 July
#
#  \author Tomas Gonzalo
#          (tomas.gonzalo@monash.edu)
#  \date 2019 Sep
#

import numpy as np
import re

from setup import *
from files import *

def get_model_parameters(parameters, partlist) :
    """
    Extracts the model (scan) parameters out of the full parameter list
    """

    model_parameters = []

    for param in parameters:
 
        # if the parameter is in MINPAR or EXTPAR or in any BLOCKIN it's not a
        # model parameter
        if param.block != "MINPAR" and param.block != "EXTPAR" :
            if not (param.block != None and param.block.endswith("IN")) : 
                model_parameters.append(param)

    return model_parameters

def get_model_par_name(paramname, parameters) :
    """
    Get the output name of the model parameter with 
    name or alt_name equal to paramname
    """
    
    for name, param in parameters.iteritems():
        if paramname == param.bcs :
            return name
    for name, param in parameters.iteritems():
        if paramname == param.name or paramname == param.alt_name :
            return param.name


def add_to_model_hierarchy(spectrum_name, model_name, model_params):
    """
    Adds a model to the model hierarchy. This means we create any
    new header files in the model directory, i.e.
    Models/include/gambit/Models/models/<new_model>.hpp
    """

    print("Writing new spectrum, {0}".format(spectrum_name))

    towrite_header = blame_gum("/// Header file for {0}".format(model_name))
    towrite_header += (
                   "#ifndef __{0}_hpp__\n"
                   "#define __{0}_hpp__\n"
                   "\n"
    ).format(model_name)
    towrite_source = blame_gum("/// Source file for {0}".format(model_name))

    module = "Models"

    towrite_header += (
    			   "#define MODEL {0}"
    			   "\n"
    			   "  START_MODEL\n"
    			   "\n"
    ).format(model_name)

    bsm_params = []

    # Don't want the SM-like Higgs mass a fundamental parameter, nor any
    # of the SM Yukawas etc, nor any Pole_Mixings.
    for p in model_params:
        if p.gb_in == "mH" and p.tag == "Pole_Mass": continue
        if p.tag == "Pole_Mixing": continue
        if p.sm == True: continue
        bsm_params.append(p)

    params = []

    for i in bsm_params:
        if i.shape == 'scalar' or i.shape == None: params.append(i.gb_in)
        elif re.match("m[2-9]x[2-9]", i.shape): 
            size = int(i.shape[-1])
            for j in xrange(size):
                for k in xrange(size):
                    params.append(i.gb_in + str(j+1) + 'x' + str(k+1))
        elif re.match("v[2-9]", i.shape): 
            size = int(i.shape[-1])
            for j in xrange(size):
                params.append(i.gb_in + str(j+1))

    # No double counting (also want to preserve the order)
    norepeats = []
    [norepeats.append(i) for i in params if not i in norepeats]


    # Chunk this up into groups of no more than 9, so the DEFINEPARS macro works.
    definepars = [norepeats[i:i+9] for i in range(0,len(norepeats),9)]

    for i in range(len(definepars)):
        towrite_header += "  DEFINEPARS({0})\n".format(', '.join(definepars[i]))

    towrite_header += (
                   "\n"
                   "#undef MODEL\n"
                   "\n"
                   "#endif\n"
    )

    return towrite_header

def find_parents_params(parent):
    """
    Returns all params in the parent model.
    """

    module = "Models"
    fname = "models/" + parent + ".hpp"
    location = full_filename(fname, module)

    lookup = "#define MODEL " + parent
    term = "#undef MODEL"

    num = find_string(fname, module, lookup)[1]

    parent_params = []

    with open(location, 'r') as f:
        for num, line in enumerate(f, 1+num):
            if "DEFINEPARS" in line:
                params = re.compile(r"\((.*)\)").search(line).group(1)
                parent_params += params.split(",")
            if term in line:
                break

    return parent_params

def find_tree_root(parent):
    """
    Traces up a model tree to find the 'root' of the model tree,
    i.e. which existing Spectrum a new model is allowed to use.
    """

    module = "Models"
    fname = "models/" + newmodel + ".hpp"

    newmodel = parent
    root = False
    lookup = "#define PARENT "

    while root == False:
        location = full_filename(fname, module)
        if lookup in open(location, 'r').read():
            lines = open(location, 'r').readlines()
            for line in lines:
              if lookup in line:
                newmodel = line.split(' ')[-1].strip('\n')
        else:
            root = True

    return newmodel

def write_spectrumcontents(gambit_model_name, model_parameters):
    """
    Writes SpectrumContents for a new model:
    Models/src/SpectrumContents/<gambit_model>.cpp.
    """

    intro = (
          "///  Class defining the parameters that SubSpectrum\n"
          "///  objects providing " + gambit_model_name + "\n"
          "///  spectrum data must provide."
    )

    towrite = blame_gum(intro)

    # TODO generate list of sizes required from the parameters.

    towrite += (
            # "#ifndef __{0}_contents_hpp__\n"
            # "#define __{0}_contents_hpp__\n\n"
            "#include \"gambit/SpecBit/RegisteredSpectra.hpp\""
            "\n"
            "\n"
            "namespace Gambit\n"
            "{{\n"
            "SpectrumContents::{0}::{0}() : Contents(\"{0}\")\n"
            "{{\n"
            "std::vector<int> scalar = initVector(1);"
            " // i.e. get(Par::Tag, \"name\")\n"
            "std::vector<int> v2     = initVector(2);"
            " // i.e. get(Par::Tag, \"name\", i)\n"
            "std::vector<int> v3     = initVector(3);   // \"\n"
            "std::vector<int> v4     = initVector(4);   // \"\n"
            "std::vector<int> v5     = initVector(5);   // \"\n"
            "std::vector<int> v6     = initVector(6);   // \"\n"
            "std::vector<int> m2x2   = initVector(2,2);"
            " // i.e. get(Par::Tag, \"name\", i, j)\n"
            "std::vector<int> m3x3   = initVector(3,3); // \"\n"
            "std::vector<int> m4x4   = initVector(4,4); // \"\n"
            "std::vector<int> m5x5   = initVector(5,5); // \"\n"
            "std::vector<int> m6x6   = initVector(6,6); // \"\n"
            "\n"
    ).format(gambit_model_name)

    # Now add each parameter to the model file.
    for i in np.arange(len(model_parameters)):

        if not isinstance(model_parameters[i], SpectrumParameter):
            raise GumError(("\n\nModel Parameters at position " + i +
                            "not passed as instance of class "
                            "SpectrumParameter."))

        mp = model_parameters[i]

        # Default shape to 'scalar'
        if mp.shape:
            shape = mp.shape
        else:
            shape = "scalar"

        # If we've extracted some information about the block, then
        # add the information to the addParameter macro
        if mp.block:
            extra = ", \"" + mp.block + "\", " + str(mp.index)
        else: extra = ""

        # Write the addParameter macro to initialise each SpectrumParameter
        # object within the SubSpectrum.

        e = mp.name[1:] if mp.tag == "Pole_Mass" else mp.name 
        towrite += (
                "addParameter(Par::{0}, \"{1}\", {2}{3});\n"
                ).format(mp.tag.replace("\"",""), 
                         e, shape, extra)


    towrite += (
            "\n"
            "} // namespace Models\n"
            "} // namespace Gambit\n"
            #"#endif\n"
    )

    contents = indent(towrite)
    return contents

def add_to_registered_spectra(gambit_model):
    """
    Adds new model entry to SpecBit/RegisteredSpectra.hpp
    """

    lookup = "Contents"
    newentry = "    struct {0:21}: Contents {{ {0}(); }};\n".format(gambit_model)
    filename = "RegisteredSpectra.hpp"
    module = "SpecBit"
    location = full_filename(filename, module)
    linenum = 0 # Position of last entry in RegisteredSpectra
    with open(location) as f:
        for num, line in enumerate(f, 1):
            if lookup in line:
                linenum = num
            if newentry in line:
                raise GumError(("\n\nModel {0} already exists in GAMBIT.").format(gambit_model))

    return newentry, linenum