"""
Master module for all Models related routines.
"""

import numpy as np
import re

from setup import *
from files import *

def add_to_model_hierarchy(new_spectrum, spectrum_name, model_name, 
                           model_params, parent=None, children=None, 
                           friend=None, 
                           translation_functions_p = None, 
                           translation_functions_c = None, 
                           translation_functions_f = None):
    """
    Adds a model to the model hierarchy. This means we create any 
    new header files in the model directory, i.e.
    Models/include/gambit/Models/models/<new_model>.hpp, and edit any
    parent/children headers. Writes translation functions etc. in
    Models/src/models/<new_model>.cpp if needed.
    """
    
    if parent or children or friend:
        new_src_file = True
    else:
        new_src_file = False
                
    if new_spectrum == True:
        print("Writing new spectrum, {0}".format(spectrum_name))
        return
      
    towrite_header = blame_gum("/// Header file for {0}".format(model_name))
    towrite_header += (
                   "#ifndef __{0}_hpp__\n"
                   "#define __{0}_hpp__\n"
                   "\n"                   
    ).format(model_name)
    towrite_source = blame_gum("/// Source file for {0}".format(model_name))
    
    module = "Models"
    header = True
    
    if parent:
        towrite_header += (
                       "#include \"gambit/Models/models/{0}.hpp\"\n"
                       "\n"
        ).format(parent)
        
    towrite_header += "#define MODEL {0}\n".format(model_name)
    
      
    if parent:
        if not find_file("models/" + parent, module, header):
            raise GumError(("\n\nParent model {0} not found. Please check " 
                            "your .gum file.").format(parent))
    
        if not translation_functions_p:
            raise GumError(("\n\nParent function specified, but no translation "
                            "function. Please check your .gum file."))
                            
        parent_params = find_parents_params(parent)
        
        towrite_header += (
                       "#define PARENT {0}\n"
                       "  START_MODEL\n"
                       "  INTERPRET_AS_PARENT_FUNCTION({1}_to_{0})\n"
        ).format(parent, model_name)
                      
        print(("No need to write a new Spectrum. The model will piggyback on "
               " parent Spectrum, {0}...").format(spectrum_name))
    else:
        print("No parent model specified. GUM will write a new Spectrum.")
    
    towrite_header += "  DEFINEPARS({0})\n".format(', '.join(model_params))
            
    if parent:
        towrite_header += "#undef PARENT\n"
    
    towrite_header += (
                   "#undef MODEL\n"
                   "\n"
                   "#endif\n"   
    )

    if new_src_file:
        towrite_source += (
                       "#include <string>\n"
                       "#include <vector>\n"
                       "\n"
                       "#include \"gambit/Models/model_macros.hpp\"\n"
                       "#include \"gambit/Models/model_helpers.hpp\"\n"
                       "#include \"gambit/Logs/logger.hpp\"\n"
                       "#include \"gambit/Utils/util_functions.hpp\"\n"
                       "\n"
                       "#include \"gambit/Models/models/{0}.hpp\"\n"
                       "#include \"gambit/Models/models/{1}.hpp\"\n"
                       "\n"
                       "using namespace Gambit::Utils;\n"
                       "\n"
                       "#define MODEL {0}\n"
                       "#define PARENT {1}\n"
                       "void MODEL_NAMESPACE::{0}_to_{1} (const "
                       "ModelParameters &myP, ModelParameters &targetP)\n"
                       "{{\n"
                       "USE_MODEL_PIPE(PARENT)\n"
                       "logger() << \"Running interpret_as_parent calculations "
                       "for {0} --> {1}...\" << LogTags::info<<EOM;\n"
                       "\n"
                       + translation_functions(model_params, parent_params,
                                               translation_functions_p) +
                       "\n"
                       "}}\n"
                       "#undef PARENT\n"
                       "#undef MODEL\n"
    ).format(model_name, parent)        
        
    #print towrite_header
    #print
    #print indent(towrite_source)
    
    return towrite_header, indent(towrite_source)

def find_parents_params(parent):
    """
    Returns all params in the parent model.
    """

    module = "Models"
    header = True
    location = full_filename("models/" + parent, module, header)
    
    lookup = "#define MODEL " + parent
    term = "#undef MODEL"
    
    num = find_string("models/" + parent, module, lookup, header)[1]
    
    parent_params = []
    
    with open(location, 'r') as f:
        for num, line in enumerate(f, 1+num):
            if "DEFINEPARS" in line:
                params = re.compile( "\((.*)\)" ).search(line).group(1)
                parent_params += params.split(",")
            if term in line:
                break
                
    return parent_params
            
 
def translation_functions(model_params, parent_params, translations):
    """
    Writes translation functions.
    """
    
    
    towrite = ""
    for i in xrange(len(model_params)):
        mp = model_params[i]
        t = ""
        if mp in translations:
            t = translations[mp]
        else:
            t = "0"
        towrite += "targetP.setValue(\"{0}\", {1});\n".format(mp, t)
    return towrite
    
    
def check_spectrum(spectrum_name, parent):
    """
    Traces up the model tree to check if there is an ancestor whose spectrum
    exists.
    """    
    
    # Hardcoded spectra. All other ones just have _spectrum attached to them.
    spectra = {'MSSM63atQ':'MSSM_spectrum'}
        
  
def find_tree_root(parent):
    """
    Traces up a model tree to find the 'root' of the model tree, 
    i.e. which existing Spectrum a new model is allowed to use.
    """

    module = "Models"
    header = True
    
    newmodel = parent
    root = False
    lookup = "#define PARENT "
    
    while root == False:
        location = full_filename("models/" + newmodel, module, header)
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
            "std::vector<int> scalar = initVector(1);"
            " // i.e. get(Par::Tag, \"name\")\n"
            "std::vector<int> v2     = initVector(2);"
            " // i.e. get(Par::Tag, \"name\", i)\n"
            "std::vector<int> v3     = initVector(3);   // \"\n"
            "std::vector<int> v4     = initVector(4);   // \"\n"
            "std::vector<int> v6     = initVector(6);   // \"\n"
            "std::vector<int> m2x2   = initVector(2,2);"
            " // i.e. get(Par::Tag, \"name\", i, j)\n"
            "std::vector<int> m3x3   = initVector(3,3); // \"\n"
            "std::vector<int> m4x4   = initVector(4,4); // \"\n"
            "std::vector<int> m6x6   = initVector(6,6); // \"\n"
            "\n"
    ).format(gambit_model_name)

    # Now add each parameter to the model file.
    for i in np.arange(len(model_parameters)):
        if not isinstance(model_parameters[i], SpectrumParameter):
            raise GumError(("\n\nModel Parameters at position " + i + 
                            "not passed as instance of class "
                            "SpectrumParameter."))
                     
    if model_parameters[i].shape:
        shape = ", " + model_parameters[i].shape
    else:
        shape = ""
    towrite += "addParameter(Par::{0}, \"{1}\"{2});\n".format(model_parameters[i].tag.replace("\"",""), model_parameters[i].name, shape)
          
    towrite += (
            "\n"
            "} // namespace Models\n"
            "} // namespace Gambit\n"
            "#endif\n"
    )
      
    contents = indent(towrite)
    return contents
      
def write_subspectrum_wrapper(gambit_model_name, model_parameters):
    """
    Writes spectrum object wrapper for new model:
    Models/include/gambit/Models/SimpleSpectra/<new_model_name>SimpleSpec.hpp.
    """
      
    # Classes make life easier 
    class SpecGetAndSet:

        def __init__(self, shape, size, param, getter, setter):
            self.shape = shape
            self.size = size
            self.param = param
            self.getter = getter
            self.setter = setter

    spectrumparameters = []

    modelSS = gambit_model_name + "SimpleSpec"
    modelclass = gambit_model_name + "Model"

    # Go through model, and create down all members of the model object,
    # all getter functions, all setter functions, all sizes...

    for i in np.arange(len(model_parameters)):
        if not isinstance(model_parameters[i], SpectrumParameter):
            raise GumError(("\n\nModel Parameters at position " + i + 
                            " not passed as instance of class "
                            "SpectrumParameter."))

        paramname = gambit_model_name + "_" + model_parameters[i].fullname.replace("\"","")

        if model_parameters[i].tag == "Pole_Mass":
            paramname += "_Pole_Mass"

        shape = "scalar"
        size = 1

        if model_parameters[i].shape:
            if re.match("v[2-9]", model_parameters[i].shape):
                shape = "vector"
                size = model_parameters[i].shape[-1]
            elif re.match("m[2-9]x[2-9]", model_parameters[i].shape):
                # Assuming all matrices will be square...
                shape = "matrix"
                size = model_parameters[i].shape[-1]

        getter = "get_" + model_parameters[i].fullname

        if model_parameters[i].tag == "Pole_Mass":
            getter += "PoleMass"

        setter = "set_" + model_parameters[i].fullname

        if model_parameters[i].tag == "Pole_Mass":
            setter += "PoleMass"

        x = SpecGetAndSet(shape, size, paramname, getter, setter)
        spectrumparameters.append(x)

    intro_message = (
            "///  A simple SubSpectrum wrapper for\n"
            "///  " + gambit_model_name + ". No RGEs included."
    )

    towrite = blame_gum(intro_message)

    towrite += (
            "#ifndef __{0}.hpp__\n"
            "#define __{0}.hpp__\n"
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
            "struct {2}\n"
            "{{\n"
    ).format(modelSS, gambit_model_name, modelclass)

    # Now add each parameter to the model file.

    for i in range(0, len(spectrumparameters)):

        sp = spectrumparameters[i]
        
        if sp.shape == "scalar":
            size = ""
        elif sp.shape == "vector":
            size = "[{0}]".format(sp.size)
        elif sp.shape == "matrix":
            size = "[{0}][{0}]".format(sp.size)

        towrite += "double {0}{1};\n".format(sp.param, size)

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
            "typedef SpectrumContents::{1}Contents;\n"
            "}};\n"
            "\n"
            "namespace Models\n"
            "{{\n"
            "class {0} : public Spec<{0}>\n"
            "{{\n"
            "private:\n {2}  params;\n"
            "typedef {0} Self;\n"
            "\n"
            "public:\n"
            "/// Constructors & destructors\n"
            "{0}(const {2}& p)\n"
            " : params(p)\n"
            "{{}}\n"
            "\n"
            "/// Wrapper functions to parameter object.\n"
    ).format(modelSS, gambit_model_name, modelclass)

    # Would be neater (here) to write get_x and set_x at the same time,
    # but following current format...

    # Getter functions
    for i in np.arange(len(spectrumparameters)):

        sp = spectrumparameters[i]
          
        if sp.shape == "scalar":
            size = ""
            indices = ""
        if sp.shape == "vector":
            size = "int i"
            indices = "[i]"
        elif sp.shape == "matrix":
            size = "int i, int j"
            indices = "[i][j]"

        towrite += "double {0}({1}) const {{ return params.{2}{3}; }}\n".format(sp.getter, size, sp.param, indices)

    towrite += "\n"

    # Setter functions
    for i in np.arange(len(spectrumparameters)):

        sp = spectrumparameters[i]
          
        if sp.shape == "scalar":
            size = ""
            indices = ""
        if sp.shape == "vector":
            size = ", int i"
            indices = "[i]"
        elif sp.shape == "matrix":
            size = ", int i, int j"
            indices = "[i][j]"

        towrite += "double {0}(double in{1}) {{ params.{2}{3}=in; }}\n".format(sp.setter, size, sp.param, indices)

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
        towrite += "typedef typename MTget::FInfo2W FInfo2W;\n"
    if v:
        towrite += "typedef typename MTget::FInfo1W FInfo1W;\n"

    # Remove all duplicates. These values tell us which indices we need to
    # include for the FInfo routines.
    sizes = list(set(sizes))
        
    for i in np.arange(len(sizes)):
        fnname = "i" + "".join(str(j) for j in np.arange(int(sizes[i])))
        
        towrite += (
                "static const {0}v[] = {{{1}}};\n"
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
            finf = " &Self::{}".format(sp.getter)
        elif sp.shape == "vector":
            size = "1"
            index = "i" + "".join(str(j) for j in np.arange(int(sp.size)))
            finf = "FInfo1W(&Self::{0}, {1})".format(sp.getter, index)
        elif sp.shape == "matrix":
            size = "2"
            index = "i" + "".join(str(j) for j in np.arange(int(sp.size)))
            finf = "FInfo2W(&Self::{0}, {1}, {1})".format(sp.getter, index) 

        towrite += (
                "getters[{0}].map{1}"
                "W[\"{2}\"] = {3};\n"
        ).format(mp.tag, size, mp.name, finf)

    towrite += (
            "\n"
            "return getters;\n"
            "}\n"
            "\n"
            "static SetterMaps fill_setter_maps()\n"
            "{\n"
            "SetterMaps setters;\n"
            "\n"
    )

    if m:
        towrite += "typedef typename MTset::FInfo2W FInfo2W;\n"
    if v:
        towrite += "typedef typename MTset::FInfo1W FInfo1W;\n"

    # Remove all duplicates. These values tell us which indices we need to
    # include for the FInfo routines.
    sizes = list(set(sizes))
    
    for i in np.arange(len(sizes)):
        fnname = "i" + "".join(str(j) for j in np.arange(int(sizes[i])))
        
        towrite += (
                "static const {0}v[] = {{{1}}};\n"
                "static const std::set<int> {0}({0}v, Utils::endA({0}v));"
                "\n"
        ).format(fnname, ",".join(str(j) for j in np.arange(int(sizes[i]))))

    towrite += "\nusing namespace Par;\n\n"

    for i in range(0, len(model_parameters)):
        sp = spectrumparameters[i]
        mp = model_parameters[i]
        
        if sp.shape == "scalar":
            size = "0"
            finf = " &Self::{}".format(sp.getter)
        elif sp.shape == "vector":
            size = "1"
            index = "i" + "".join(str(j) for j in np.arange(int(sp.size)))
            finf = "FInfo1W(&Self::{0}, {1})".format(sp.setter, index)
        elif sp.shape == "matrix":
            size = "2"
            index = "i" + "".join(str(j) for j in np.arange(int(sp.size)))
            finf = "FInfo2W(&Self::{0}, {1}, {1})".format(sp.setter, index) 
        
        towrite += (
                "setters[{0}].map{1}"
                "W[\"{2}\"] = {3};\n"
        ).format(mp.tag, size, mp.name, finf)

    towrite += (
            "\n"
            "return setters;\n"
            "}\n"
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
    newentry = "    struct {0:16}: SubSpectrumContents {{ {0}(); }};\n".format(gambit_model)
    filename = "SpectrumContents/RegisteredSpectra"
    module = "Models"
    location = full_filename(filename, module, header=True)
    linenum = 0 # Position of last entry in RegisteredSpectra
    with open(location) as f:
        for num, line in enumerate(f, 1):
            if lookup in line:
                linenum = num
            if newentry in line:
                raise GumError(("\n\nModel {0} already exists in GAMBIT.").format(gambit_model))
    
    return newentry, linenum          
    amend_file(filename, module, newentry, linenum, header=True)
