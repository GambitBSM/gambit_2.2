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

def get_model_parameters(parameters, partlist):
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


    # Replace all trailing + and - with pm
    for param in parameters:
        param.name = re.sub(r'(.*)[-+]', r'\1pm', param.name)

    return model_parameters

def get_model_par_name(paramname, parameters):
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


    towrite += (
            "#ifndef __{0}_contents_hpp__\n"
            "#define __{0}_contents_hpp__\n\n"
            "#include \"gambit/Models/SpectrumContents/RegisteredSpectra.hpp\""
            "\n"
            "\n"
            "namespace Gambit\n"
            "{{\n"
            "SpectrumContents::{0}::{0}()\n"
            "{{\n"
            "setName(\"{0}\");\n"
            "\n"
    ).format(gambit_model_name)

    # Extract all shapes from here
    shapes = sorted(list(set([mp.shape for mp in model_parameters])), 
                             reverse=True)

    # Scalars first...
    if "scalar" in shapes:
        towrite += (
                "std::vector<int> scalar = initVector(1);"
                " // i.e. get(Par::Tag, \"name\")\n"
        )

    for shape in shapes:
        if re.search(r'v\d', shape):
            towrite += (
                    "std::vector<int> {0}     = initVector({1});"
                    " // i.e. get(Par::Tag, \"name\", i)\n"
            ).format(shape, shape[1:])
            continue
        elif re.search(r'm\dx\d', shape):
            towrite += (
                    "std::vector<int> {0}  = initVector({1},{2});"
                    " // i.e. get(Par::Tag, \"name\", i, j)\n"
            ).format(shape, shape[1], shape[3])
            continue
        elif shape == "scalar": continue
        else:
            raise GumError(("Shape {0} not recognised...").format(shape))

    towrite += "\n"

    # Now add each parameter to the model file.
    for i in np.arange(len(model_parameters)):

        mp = model_parameters[i]
        
        if not isinstance(mp, SpectrumParameter):
            raise GumError(("\n\nModel Parameters at position " + i +
                            "not passed as instance of class "
                            "SpectrumParameter."))

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

        e = mp.fullparticlename if mp.tag == "Pole_Mass" else mp.name
        towrite += (
                "addParameter(Par::{0}, \"{1}\", {2}{3});\n"
                ).format(mp.tag.replace("\"",""), 
                         e, shape, extra)


    towrite += (
            "\n"
            "} // namespace Models\n"
            "} // namespace Gambit\n"
            "#endif\n"
    )

    contents = indent(towrite)
    return contents

def write_subspectrum_wrapper(gambit_model_name, model_parameters,
                              bsm_partlist, mixings, gambit_pdgs):
    """
    Writes spectrum object wrapper for new model:
    Models/include/gambit/Models/SimpleSpectra/<new_model_name>SimpleSpec.hpp.
    """

    # Classes make life easier
    class SpecGetAndSet:

        def __init__(self, shape, size, param, getter, setter, block, index):
            self.shape = shape
            self.size = size
            self.param = param
            self.getter = getter
            self.setter = setter
            self.block = block
            self.index = index

    spectrumparameters = []

    modelSS = gambit_model_name + "SimpleSpec"
    modelclass = gambit_model_name + "Model"

    # Go through model, and create down all members of the model object,
    # all getter functions, all setter functions, all sizes...

    # Keep track of which masses we've added
    addedpdgs = []

    for i in np.arange(len(model_parameters)):
        
        par = model_parameters[i]
        
        if not isinstance(par, SpectrumParameter):
            raise GumError(("\n\nModel Parameters at position " + i +
                            " not passed as instance of class "
                            "SpectrumParameter."))

        if par.sm:
            e = ""
        else:
            e = gambit_model_name + "_"

        # Remove the trailing 'm' for a Pole_Mass
        if par.tag == "Pole_Mass":
            paramname = e + par.fullname[1:].strip('~') + "_Pole_Mass"
        else:
            paramname = e + par.fullname

        shape = "scalar"
        size = 1

        if par.shape:
            if re.match("v[2-9]", par.shape):
                shape = "vector"
                size = par.shape[-1]
            elif re.match("m[2-9]x[2-9]", par.shape):
                # Assuming all matrices will be square...
                shape = "matrix"
                size = par.shape[-1]

        if par.tag == "Pole_Mass":
            setter = "set_" + par.fullname[1:].strip('~') + "PoleMass"
            getter = "get_" + par.fullname[1:].strip('~') + "PoleMass"
        else:
            setter = "set_" + par.fullname
            getter = "get_" + par.fullname

        # Replace all plusses and minuses with 'pm'
        setter = setter.replace("-","pm").replace("+","pm")
        getter = getter.replace("-","pm").replace("+","pm")
        paramname = paramname.replace("-","pm").replace("+","pm")

        # Get the block and index if relevant
        block = par.block
        index = par.index

        # Save PDG codes so we don't double count
        if par.tag == "Pole_Mass":
            addedpdgs.append(par.index)

        x = SpecGetAndSet(shape, size, paramname, getter, setter, block, index)
        spectrumparameters.append(x)

    # Go through BSM particle list and add the pole masses to the list of 
    # params a spectrum object should interface to 
    for particle in bsm_partlist:

        # Don't double count particles!
        if particle.PDG_code in addedpdgs: 
            continue

        pname = particle.name.strip('~')

        paramname = gambit_model_name + "_" + pname + "_Pole_Mass"
        setter = "set_" + pname + "PoleMass"
        getter = "get_" + pname + "PoleMass"

        # Replace all plusses and minuses with 'pm'
        setter = setter.replace("-","pm").replace("+","pm")
        getter = getter.replace("-","pm").replace("+","pm")
        paramname = paramname.replace("-","pm").replace("+","pm")

        block = "MASS"
        index = abs(particle.PDG_code)

        x = SpecGetAndSet("scalar", 1, paramname, getter, setter, block, index)
        spectrumparameters.append(x)

    # Convert the C++ dict to python properly
    mixingdict = dict((m.key(),m.data()) for m in mixings)

    # print mixingdict
    # TODO mixing dict
    # # Same for the mixings
    # for mix in mixingdict:
    #     ...

    intro_message = (
            "///  A simple SubSpectrum wrapper for\n"
            "///  " + gambit_model_name + ". No RGEs included."
    )

    towrite = blame_gum(intro_message)

    towrite += (
            "#ifndef __{0}_hpp__\n"
            "#define __{0}_hpp__\n"
            "\n"
            "#include \"gambit/Elements/spec.hpp\"\n"
            "#include \"gambit/Models/SpectrumContents/"
            "RegisteredSpectra.hpp\"\n"
            "\n"
            "namespace Gambit\n"
            "{{\n"
            "namespace Models\n"
            "{{\n"
            "/// Simple {1} model object.\n"
            "class {2} : public SLHAeaModel\n"
            "{{\n"
            "\n"
            "public:\n"
            "  /// @{{ Constructors\n"
            "{2}(const SLHAstruct &input)\n"
            " : SLHAeaModel(input)\n"
            "{{}}\n"
            "  /// @}}\n"
            "\n"
            "  /// @{{ Getters for {1} information\n"
    ).format(modelSS, gambit_model_name, modelclass)

    # Now add each parameter to the model file.
    # TODO: TG: In the MSSM-like SimpleSpec stuff, the getters go here, not the parameters, so delete if appropriate

    #for i in range(0, len(spectrumparameters)):

    #    sp = spectrumparameters[i]

    #    if sp.shape == "scalar":
    #        size = ""
    #    elif sp.shape == "vector":
    #        size = "[{0}]".format(sp.size)
    #    elif sp.shape == "matrix":
    #        size = "[{0}][{0}]".format(sp.size)

    #    towrite += "double {0}{1};\n".format(sp.param, size)

    # Getter functions
    for i in np.arange(len(spectrumparameters)):

        sp = spectrumparameters[i]

        if sp.shape == "scalar":
            size = ""
            indices = sp.index
        if sp.shape == "vector":
            size = "int i"
            indices = "i"
        elif sp.shape == "matrix":
            size = "int i, int j"
            indices = "i,j"

        towrite += (
            "double {0}({1}) const {{ return getdata(\"{2}\",{3}); }}\n"
        ).format(sp.getter, size, sp.block, indices)

    towrite += "  /// @}}\n\n"

    towrite += (
            "}};\n"
            "\n"
            "/// Forward declare the wrapper class so that we can use it\n"
            "/// as the template parameter for the SpecTraits specialisation.\n"
            "class {0};\n"
            "}}"
            "\n"
            "\n"
            "/// Specialisation of traits class needed to inform "
            "base spectrum class of the Model and Input types\n"
            "template <> \n"
            "struct SpecTraits<Models::{0}> : DefaultTraits\n"
            "{{\n"
            "static std::string name() {{ return \"{0}\"; }}\n"
            "typedef SpectrumContents::{1} Contents;\n"
            "typedef Models::{2} Model;\n"
            "}};\n"
            "\n"
            "namespace Models\n"
            "{{\n"
            "class {0} : public SLHASimpleSpec<{0}>\n"
            "{{\n"
            "\n"
            "public:\n"
            "  /// @{{\n"
            "/// Constructor via SLHAea object\n"
            "{0}(const SLHAea::Coll& input)\n"
            " : SLHASimpleSpec(input)\n"
            "{{}}\n"
            "\n"
            "/// Copy constructor\n"
            "{0}(const {0}& other)\n"
            " : SLHASimpleSpec(other)\n"
            "{{}}\n"
            "\n"
            "/// Destructor\n"
            "virtual ~{0}() {{}};\n"
            "\n"
            "static int index_offset() {{return -1;}}\n"
            "\n"
            "/// Construct the SubSpectrumContents\n"
            "const SpectrumContents::{1} contents;\n"
            "\n"
            "/// Add SLHAea object using the SimpleSpec_to_SLHAea routine\n"
            "void add_to_SLHAea(int /*slha_version*/, SLHAea::Coll& slha) const\n"
            "{{\n"
            "// Add SPINFO data if not already present\n"
            "SLHAea_add_GAMBIT_SPINFO(slha);\n"
            "\n"
            "// All blocks given in the SimpleSpec\n"
            "\nadd_SimpleSpec_to_SLHAea(*this, slha, contents);\n"
            "}}\n"
            "\n" 
            "/// Wrapper functions to parameter object.\n"
    ).format(modelSS, gambit_model_name, modelclass)

    # Would be neater (here) to write get_x and set_x at the same time,
    # but following current format...
    # TODO: TG: In the MSSM-like SimpleSpec stuff, getters go on the model.
    # Delete this if confirmed

    # Getter functions
    #for i in np.arange(len(spectrumparameters)):

    #    sp = spectrumparameters[i]

    #    if sp.shape == "scalar":
    #        size = ""
    #        indices = ""
    #    if sp.shape == "vector":
    #        size = "int i"
    #        indices = "[i]"
    #    elif sp.shape == "matrix":
    #        size = "int i, int j"
    #        indices = "[i][j]"

    #    towrite += "double {0}({1}) const {{ return params.{2}{3}; }}\n".format(sp.getter, size, sp.param, indices)

    #towrite += "\n"

    # TODO: TG: I don't think setters are needed, but check. If so add them where the model definition is
    # Setter functions
    #for i in np.arange(len(spectrumparameters)):

    #    sp = spectrumparameters[i]

    #    if sp.shape == "scalar":
    #        size = ""
    #        indices = ""
    #    if sp.shape == "vector":
    #        size = ", int i"
    #        indices = "[i]"
    #    elif sp.shape == "matrix":
    #        size = ", int i, int j"
    #        indices = "[i][j]"

    #    towrite += "void {0}(double in{1}) {{ params.{2}{3}=in; }}\n".format(sp.setter, size, sp.param, indices)

    towrite += (
            "\n"
            "/// Map fillers\n"
            "static GetterMaps fill_getter_maps()\n"
            "{\n"
            "GetterMaps getters;\n"
            "\n"
    )

    # Add necessary function pointer maps.
    v = False
    m = False
    sizes = []
    sizes = []
    for i in range(0, len(spectrumparameters)):
        sp = spectrumparameters[i]
        if sp.shape == "vector":
            v = True
            sizes.append(sp.size)
        elif sp.shape == "matrix":
            m = True
            sizes.append(sp.size)

    if m:
        towrite += "typedef typename MTget::FInfo2 FInfo2;\n"
    if v:
        towrite += "typedef typename MTget::FInfo1 FInfo1;\n"

    # Remove all duplicates. These values tell us which indices we need to
    # include for the FInfo routines.
    sizes = list(set(sizes))

    for i in np.arange(len(sizes)):
        fnname = "i" + "".join(str(j) for j in np.arange(int(sizes[i])))

        towrite += (
                "static const int {0}v[] = {{{1}}};\n"
                "static const std::set<int> {0}({0}v, Utils::endA({0}v));"
                "\n"
        ).format(fnname, ",".join(str(j) for j in np.arange(int(sizes[i]))))

    towrite += "\nusing namespace Par;\n\n"

    # Now add getter maps
    for i in np.arange(len(model_parameters)):
        sp = spectrumparameters[i]
        mp = model_parameters[i]

        if sp.shape == "scalar":
            size = "0"
            finf = " &Model::{}".format(sp.getter)
        elif sp.shape == "vector":
            size = "1"
            index = "i" + "".join(str(j) for j in np.arange(int(sp.size)))
            finf = "FInfo1(&Model::{0}, {1})".format(sp.getter, index)
        elif sp.shape == "matrix":
            size = "2"
            index = "i" + "".join(str(j) for j in np.arange(int(sp.size)))
            finf = "FInfo2(&Model::{0}, {1}, {1})".format(sp.getter, index)

        e = mp.fullparticlename if mp.tag == "Pole_Mass" else mp.name

        towrite += (
                "getters[{0}].map{1}"
                "[\"{2}\"] = {3};\n"
        ).format(mp.tag, size, e, finf)

    for particle in bsm_partlist:

        # Don't double count particles!
        if particle.PDG_code in addedpdgs: 
            continue

        towrite += (
                "getters[\"Pole_Mass\"].map0[\"{0}\"] = "
                "&Model::get_{1}PoleMass;\n"
        ).format(pdg_to_particle(particle.PDG_code, gambit_pdgs), 
                 particle.name.strip('~'))


    towrite += (
            "\n"
            "return getters;\n"
            "}\n"
            "\n"
    )
    # TODO: TG: I'm not sure setter maps are needed, uncomment if so
    #        "static SetterMaps fill_setter_maps()\n"
    #        "{\n"
    #        "SetterMaps setters;\n"
    #        "\n"
    #)

    #if m:
    #    towrite += "typedef typename MTset::FInfo2W FInfo2W;\n"
    #if v:
    #    towrite += "typedef typename MTset::FInfo1W FInfo1W;\n"

    # Remove all duplicates. These values tell us which indices we need to
    # include for the FInfo routines.
    #sizes = list(set(sizes))

    #for i in np.arange(len(sizes)):
    #    fnname = "i" + "".join(str(j) for j in np.arange(int(sizes[i])))

    #    towrite += (
    #            "static const int {0}v[] = {{{1}}};\n"
    #            "static const std::set<int> {0}({0}v, Utils::endA({0}v));"
    #            "\n"
    #    ).format(fnname, ",".join(str(j) for j in np.arange(int(sizes[i]))))

    #towrite += "\nusing namespace Par;\n\n"

    #for i in range(0, len(model_parameters)):
    #    sp = spectrumparameters[i]
    #    mp = model_parameters[i]

    #    if sp.shape == "scalar":
    #        size = "0"
    #        finf = " &Self::{}".format(sp.setter)
    #    elif sp.shape == "vector":
    #        size = "1"
    #        index = "i" + "".join(str(j) for j in np.arange(int(sp.size)))
    #        finf = "FInfo1W(&Self::{0}, {1})".format(sp.setter, index)
    #    elif sp.shape == "matrix":
    #        size = "2"
    #        index = "i" + "".join(str(j) for j in np.arange(int(sp.size)))
    #        finf = "FInfo2W(&Self::{0}, {1}, {1})".format(sp.setter, index)

    #    e = mp.fullparticlename if mp.tag == "Pole_Mass" else mp.name

    #    towrite += (
    #            "setters[{0}].map{1}"
    #            "W[\"{2}\"] = {3};\n"
    #    ).format(mp.tag, size, e, finf)

    #towrite += (
    #        "\n"
    #        "return setters;\n"
    #        "}\n"
    towrite += (
            "};\n"
            "}\n"
            "} // namespace Gambit\n"
            "#endif\n"
    )

    contents = indent(towrite)
    return contents

def add_to_registered_spectra(gambit_model):
    """
    Adds new model entry to RegisteredSpectra.hpp
    """

    lookup = "SubSpectrumContents"
    newentry = "    struct {0:21}: SubSpectrumContents {{ {0}(); }};\n".format(gambit_model)
    filename = "SpectrumContents/RegisteredSpectra.hpp"
    module = "Models"
    location = full_filename(filename, module)
    linenum = 0 # Position of last entry in RegisteredSpectra
    with open(location) as f:
        for num, line in enumerate(f, 1):
            if lookup in line:
                linenum = num
            if newentry in line:
                raise GumError(("\n\nModel {0} already exists in GAMBIT.").format(gambit_model))

    return newentry, linenum
