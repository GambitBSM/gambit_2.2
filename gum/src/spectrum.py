"""
Master module for all SpecBit related routines.
"""

from setup import *
from files import *

def write_basic_spectrum(gambit_model_name, model_parameters, spec,
                         simple_SMinputs=False, FS=False):
    """
    Writes basic spectrum object wrapper for new model:
    SpecBit/src/SpecBit_<new_model_name>. 
    """
    
    modelSS = gambit_model_name + "SimpleSpec"
    modelclass = gambit_model_name + "Model"
    modelcont = gambit_model_name + "model"
    
    intro_message = (
            "///  Implementation of SpecBit routines for \n"
            "///  " + gambit_model_name + "."
    )

    towrite = blame_gum(intro_message)

    towrite += (
            "#include \"gambit/Elements/gambit_module_headers.hpp\"\n"
            "#include \"gambit/Elements/spectrum.hpp\"\n"
            "#include \"gambit/Utils/stream_overloads.hpp\"\n"
            "#include \"gambit/Utils/util_macros.hpp\"\n"
            "#include \"gambit/SpecBit/SpecBit_rollcall.hpp\"\n"
            "#include \"gambit/SpecBit/SpecBit_helpers.hpp\"\n"
            "#include \"gambit/SpecBit/QedQcdWrapper.hpp\"\n"
            "#include \"gambit/Models/SimpleSpectra/{0}.hpp\"\n"
            "#include \"gambit/SpecBit/{1}Spec.hpp\"\n"
    ).format(modelSS, gambit_model_name)
    
    # Add FlexibleSUSY headers if flagged
    # TODO - wait until template BOSS...
    if FS:
      towrite += (
            "#include \"gambit/SpecBit/model_files_and_boxes.hpp\"\n"
      )
          
    towrite += (
            "\n"
            "namespace Gambit\n"
            "{{\n"
            "\n"
            "namespace SpecBit\n"
            "{{\n"
            "using namespace LogTags;\n"
            "\n"
            "/// Get a simple wrapper for Spectrum object.\n"
            "void get_{0}(Spectrum& result)\n"
            "{{\n"
            "namespace myPipe = Pipes::get_{0};\n"
            "const SMInputs& sminputs = *myPipe::Dep::SMINPUTS;\n"
            "\n"
            "// Initialise model object \n"
            "Models::{1} {2};\n\n"
            "// BSM parameters\n"
    ).format(spec, modelclass, modelcont)
    
    
    # Now add each BSM model parameter to spectrum
    for i in range(0, len(model_parameters)):
    
        par = model_parameters[i] 
      
        if not isinstance(par, SpectrumParameter):
            raise GumError(("\n\nModel Parameters at position " + i + 
                            "not passed as instance of class "
                            "SpectrumParameter."))
          
        if not par.sm:
            toadd = ""
            if par.tag == "Pole_Mass":
                toadd = "_Pole_Mass"
            towrite += "{0}.{1}_{2}{3} = *myPipe::Param.at(\"{4}\");\n".format(modelcont, gambit_model_name, par.fullname, toadd, par.gb_in)
                        
    if simple_SMinputs:
        towrite += add_simple_sminputs(modelcont)
      
    towrite += (
            "\n"
            "// Create a SubSpectrum object wrapper\n"
            "Models::{0} spec({1});\n\n" 
            "// Retrieve any mass cuts\n"
            "static const Spectrum::mc_info mass_cut = "
            "myPipe::runOptions->getValueOrDef<Spectrum::mc_info>"
            "(Spectrum::mc_info(), \"mass_cut\");\n"
            "static const Spectrum::mr_info mass_ratio_cut = "
            "myPipe::runOptions->getValueOrDef<Spectrum::mr_info>"
            "(Spectrum::mr_info(), \"mass_ratio_cut\");\n\n"
            "// We don't supply a LE subspectrum here; an "
            "SMSimpleSpec will therefore be automatically created "
            "from 'sminputs'\n"
            "result = Spectrum(spec,sminputs,&myPipe::Param,"
            "mass_cut,mass_ratio_cut);\n"
            "}}\n\n"
            "void get_{2}_spectrum_as_map" 
            "(std::map<std::string,double>& specmap)\n"
            "{{\n"
            "namespace myPipe = Pipes::get_{3}_as_map;\n"
            "const Spectrum& spec(*myPipe::Dep::{3};\n"
            "fill_map_from_{2}spectrum(specmap, spec);\n"
            "}}\n\n"
            "void fill_map_from_{2}spectrum"
            "(std::map<std::string, double>& specmap, "
            "const Spectrum& spec)\n"
            "{{\n"
            "/// Use SpectrumContents routines to automate\n"
            "static const SpectrumContents::{2} contents;\n"
            "static const std::vector<SpectrumParameter> "
            "required_parameters = contents.all_parameters();\n\n"
            "for(std::vector<SpectrumParameter>::const_iterator "
            "it = required_parameters.begin(); "
            "it != required_parameters.end(); ++it)\n"
            "{{\n"
            "const Par::Tags        tag   = it->tag();\n"
            "const std::string      name  = it->name();\n"
            "const std::vector<int> shape = it->shape();\n"
            "\n"
            "// Scalar case\n"
            "if(shape.size()==1 and shape[0]==1)\n"
            "{{\n"
            "std::ostringstream label;\n"
            "label << name <<" "<< Par::toString.at(tag);\n"
            "specmap[label.str()] = spec.get_HE().get(tag,name);\n"
            "}}\n"
            "// Vector case\n"
            "else if(shape.size()==1 and shape[0]>1)\n"
            "{{\n"
            "for(int i = 1; i<=shape[0]; ++i)\n"
            "{{\n"
            "std::ostringstream label;\n"
            "label << name <<\"_\"<<i<<\" \"<< Par::toString.at(tag);\n"
            "specmap[label.str()] = spec.get_HE().get(tag,name,i);\n"
            "}}\n"
            "}}\n"
            "// Matrix case\n"
            "else if(shape.size()==2)\n"
            "{{\n"
            "for(int i = 1; i<=shape[0]; ++i)\n"
            "{{\n"
            "for(int j = 1; j<=shape[0]; ++j)\n"
            "{{\n"
            "std::ostringstream label;\n"
            "label << name <<\"_(\"<<i<<\",\"<<j<<\") \""
            "<<Par::toString.at(tag);\n"
            "specmap[label.str()] = spec.get_HE().get(tag,name,i,j);\n"
            "}}\n"
            "}}\n"
            "}}\n"
            "// Deal with all other cases\n"
            "else\n"
            "{{\n"
            "// ERROR\n"
            "std::ostringstream errmsg;\n"
            "errmsg << \"Invalid parameter received while converting " 
            "{3} to map of strings!\"\n"
            "errmsg << \"Problematic parameter was: \"<< tag "
            "<<\", \" << name << \", shape=\"<< shape;\n"
            "utils_error().forced_throw(LOCAL_INFO,errmsg.str());\n"
            "}}\n"
            "}}\n"
            "}}\n\n"
            "}}\n\n"
            "}}\n"
    ).format(modelSS, modelcont, gambit_model_name, spec)
        
    filename = "SpecBit_" + gambit_model_name
    module = "SpecBit"
    contents = indent(towrite)
    
    return contents
    
    
def write_spectrum_header(model_name):
    """
    Writes the header for spectrum object,
    SpecBit/include/gambit/SpecBit/SpecBit_<model>_rollcall.hpp
    """
    
    towrite = blame_gum(("///  Rollcall declarations for routines declared \n"
                         "///  in SpecBit_{0}.cpp.".format(model_name)))
    
    towrite += (
            "#idndef __SpecBit_{0}_hpp\n"
            "#define __SpecBit_{0}_hpp\n"
            "\n"
            "  // Spectrum object\n"
            "  #define CAPABILITY {0}_spectrum\n"
            "  START_CAPABILITY\n"
            "\n"
            "    // Create simple object from SMInputs & new params.\n"
            "    #define FUNCTION get_{0}_spectrum\n"
            "    START_FUNCTION(Spectrum)\n"
            "    DEPENDENCY(SMINPUTS, SMInputs)\n"
            "    ALLOW_MODELS({0})\n"
            "    #undef FUNCTION\n"
            "\n"
            "    // Map for Spectrum, for printing.\n"
            "    #define FUNCTION get_{0}_spectrum_as_map\n"
            "    START_FUNCTION(map_str_dbl)\n"
            "    DEPENDENCY({0}_spectrum, Spectrum)\n"
            "    #undef FUNCTION\n"
            "\n"
            "  #undef CAPABILITY\n"
            "\n"
            "#endif\n"
    ).format(model_name)
    
    return towrite
    
def write_specbit_rollcall(model_name):
    """
    Writes the entry for the SpecBit rollcall header.
    """
    
    towrite = (
            "\n"
            "/// Module function declarations for SpecBit_{0}.cpp\n"
            "#include \"gambit/SpecBit/SpecBit_{0}_rollcall.hpp\"\n"
            "\n"
    ).format(model_name)
    
    return dumb_indent(2, towrite)
    
def add_simple_sminputs(model):
    """
    Adds simple SMInputs definitions to a spectrum object.
    """
      
    towrite = (
            "\n"
            "// quantities needed to fill container spectrum\n"
            "double alpha_em = 1.0 / sminputs.alphainv;\n"
            "double C = alpha_em * Pi / (sminputs.GF * pow(2,0.5));\n"
            "double sinW2 = 0.5 - pow( 0.25 - C/pow(sminputs.mZ,2) , 0.5);\n"
            "double cosW2 = 0.5 + pow( 0.25 - C/pow(sminputs.mZ,2) , 0.5);\n"
            "double e = pow( 4*Pi*( alpha_em ),0.5);\n"
            "double vev = 1. / sqrt(sqrt(2.)*sminputs.GF);\n"
            "\n"
            "// Gauge couplings\n"
            "{0}.HiggsVEV = vev;\n"    
            "{0}.g1 = e / sqrt(sinW2);\n"    
            "{0}.g2 = e / sqrt(cosW2);\n"    
            "{0}.g3 = pow( 4*Pi*( sminputs.alphaS ),0.5);\n"    
            "\n"
            "// Yukawas\n"
            "double sqrt2v = pow(2.0,0.5)/vev;\n"
            "{0}.Yu[0] = sqrt2v * sminputs.mU;\n"
            "{0}.Yu[1] = sqrt2v * sminputs.mCmC;\n"
            "{0}.Yu[2] = sqrt2v * sminputs.mT;\n"
            "{0}.Ye[0] = sqrt2v * sminputs.mE;\n"
            "{0}.Ye[1] = sqrt2v * sminputs.mMu;\n"
            "{0}.Ye[2] = sqrt2v * sminputs.mTau;\n"
            "{0}.Yd[0] = sqrt2v * sminputs.mD;\n"
            "{0}.Yd[1] = sqrt2v * sminputs.mS;\n"
            "{0}.Yd[2] = sqrt2v * sminputs.mBmB;\n"
    ).format(model)
      
    return towrite
    
