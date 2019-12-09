"""
Master module for all SpecBit related routines.
"""

from setup import *
from files import *
from cmake_variables import *

def write_spectrum(gambit_model_name, model_parameters, spec,
                   add_higgs, with_spheno, gambit_pdgs,
                   neutral_higgses, charged_higgses, blockdict, 
                   particles):
    """
    Writes the spectrum object wrapper for new model:
    SpecBit/src/SpecBit_<new_model_name>.cpp.

    GUM will create a basic spectrum container unless the user has
    specified SPheno output. If so the spectrum object will be built by 
    interfacing to SPheno. This also contains routines to interface to
    HiggsBounds/HiggsSignals via the GAMBIT HiggsCouplingTable.
    """

    modelSS = gambit_model_name + "SimpleSpec"
    modelclass = gambit_model_name + "Model"
    modelcont = gambit_model_name + "model"

    is_susy = True if "MSSM" in gambit_model_name else False
        
    intro_message = (
            "///  Implementation of SpecBit routines for \n"
            "///  " + gambit_model_name + "."
    )

    towrite = blame_gum(intro_message)

    # Simple spectrum wrapper for a new model with no SPheno interface.
    if not with_spheno:

        towrite += (
                "#include <string>\n"
                "#include <sstream>\n"
                "\n"
                "#include \"gambit/Elements/gambit_module_headers.hpp\"\n"
                "#include \"gambit/Elements/spectrum.hpp\"\n"
                "#include \"gambit/Utils/stream_overloads.hpp\"\n"
                "#include \"gambit/Utils/util_macros.hpp\"\n"
                "\n"
                "#include \"gambit/SpecBit/SpecBit_rollcall.hpp\"\n"
                "#include \"gambit/SpecBit/SpecBit_helpers.hpp\"\n"
                "#include \"gambit/Models/SpectrumContents/RegisteredSpectra.hpp\""
                "\n"
                "#include \"gambit/SpecBit/QedQcdWrapper.hpp\"\n"
                "#include \"gambit/Models/SimpleSpectra/{0}.hpp\"\n"
        ).format(modelSS, gambit_model_name)

        # If there's a Higgs dependency 
        if add_higgs:
            towrite += (
                    "#include \"gambit/Models/SimpleSpectra/"
                    "SMHiggsSimpleSpec.hpp\"\n"
            )
        
        towrite += (
                "namespace Gambit\n"
                "{{\n"
                "\n"
                "namespace SpecBit\n"
                "{{\n"
                "using namespace LogTags;\n"
                "\n"
                "/// Get a (simple) Spectrum object wrapper for {0} model.\n"
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

            # Don't have anything that's an output of spectrum computation as 
            # a scan parameter
            if par.is_output: continue
          
            if not isinstance(par, SpectrumParameter):
                raise GumError(("\n\nModel Parameters at position " + i + 
                                "not passed as instance of class "
                                "SpectrumParameter."))
              
            if not par.sm:
                
                # If it's a pole mass then append this information onto the parameter name
                toadd = "_Pole_Mass" if par.tag == "Pole_Mass" else ""
                shape = "scalar"
                size = 1

                if model_parameters[i].shape:
                    if re.match("v[2-9]", par.shape):
                        shape = "vector"
                        size = par.shape[-1]
                    elif re.match("m[2-9]x[2-9]", par.shape):
                        # Assuming all matrices will be square...
                        shape = "matrix"
                        size = par.shape[-1]

                e = par.fullname[1:].strip('~') if par.tag == "Pole_Mass" else par.fullname

                # If it's a scalar shape, just add it one by one
                if shape == "scalar":
                    towrite += (
                            "{0}.{1}_{2}{3} = *myPipe::Param[\"{4}\"];\n"
                    ).format(modelcont, gambit_model_name, e, toadd, par.gb_in)

                # If it's a matrix then do each element individually
                elif shape == "matrix":
                    for i in xrange(int(size)):
                        for j in xrange(int(size)):
                            towrite += (
                                "{0}.{1}_{2}{3}[{4}][{5}] = *myPipe::Param[\"{6}{7}x{8}\"];\n"
                            ).format(modelcont, gambit_model_name, e, toadd, i, j, par.gb_in, i+1, j+1)

                # Same deal for a vector
                elif shape == "vector":
                  for i in xrange(int(size)):
                      towrite += (
                              "{0}.{1}_{2}{3}[{4}] = *myPipe::Param[\"{5}{6}\"];\n"
                          ).format(modelcont, gambit_model_name, e, toadd, i, par.gb_in, i+1)

                else:
                    raise GumError("Parameter with shape " + shape + " is not currently supported.")

        if add_higgs:
            towrite += add_simple_sminputs(modelcont)

        towrite += (
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
        ).format(modelSS, modelcont)

    # If we have SPheno, make the spectrum use it.
    # TODO SpecBit/include/NEWMODELSpec.hpp does not get generated by GUM (yet)
    else:
        towrite += (
                "#include <cmath>\n"
                "#include \"gambit/Elements/gambit_module_headers.hpp\"\n"
                "#include \"gambit/Elements/spectrum_factories.hpp\"\n"
                "#include \"gambit/Elements/smlike_higgs.hpp\"\n"
                "#include \"gambit/Models/SimpleSpectra/{0}SimpleSpec.hpp\"\n"
                "#include \"gambit/SpecBit/{0}Spec.hpp\"\n"
                "#include \"gambit/SpecBit/SpecBit_rollcall.hpp\"\n"
                "#include \"gambit/SpecBit/SpecBit_helpers.hpp\"\n"
                "\n"
                "namespace Gambit\n"
                "{{\n"
                "\n"
                "namespace SpecBit\n"
                "{{\n"
                "using namespace LogTags;\n"
                "\n"
        ).format(gambit_model_name)

        towrite += write_spheno_spectrum_src(gambit_model_name, is_susy)

        # If we have a non-SM like number of Higgses then write an interface
        # to HiggsCouplingTable via SPheno 
        if len(neutral_higgses+charged_higgses) > 1:
            towrite += write_spheno_higgsbounds_interface(gambit_model_name, 
                                                          gambit_pdgs,
                                                          neutral_higgses,
                                                          charged_higgses)

    # Printing, fill_map_from_spectrum, and wrap it up. 
    towrite += (
            "// Declaration: print spectrum out\n"
            "void fill_map_from_{2}_spectrum(std::map<std::string,double>&, "
            "const Spectrum&);\n\n"
            "void get_{2}_spectrum_as_map" 
            "(std::map<std::string,double>& specmap)\n"
            "{{\n"
            "namespace myPipe = Pipes::get_{3}_as_map;\n"
            "const Spectrum& spec(*myPipe::Dep::{3});\n"
            "fill_map_from_{2}_spectrum(specmap, spec);\n"
            "}}\n\n"
            "void fill_map_from_{2}_spectrum"
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
            "label << name <<\" \"<< Par::toString.at(tag);\n"
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
            "{3} to map of strings!\";\n"
            "errmsg << \"Problematic parameter was: \"<< tag "
            "<<\", \" << name << \", shape=\"<< shape;\n"
            "utils_error().forced_throw(LOCAL_INFO,errmsg.str());\n"
            "}}\n"
            "}}\n"
            "}}\n\n"
            "}}\n\n"
            "}}\n"
    ).format(modelSS, modelcont, gambit_model_name, spec)
            
    return indent(towrite)
    
    
def write_spectrum_header(model_name, add_higgs, with_spheno, higgses):
    """
    Writes the header for spectrum object,
    SpecBit/include/gambit/SpecBit/SpecBit_<model>_rollcall.hpp
    """
    
    towrite = blame_gum(("///  Rollcall declarations for routines declared \n"
                         "///  in SpecBit_{0}.cpp.".format(model_name)))
    
    towrite += (
            "#ifndef __SpecBit_{0}_hpp__\n"
            "#define __SpecBit_{0}_hpp__\n"
            "\n"
            "  // Spectrum object\n"
            "  #define CAPABILITY {0}_spectrum\n"
            "  START_CAPABILITY\n"
            "\n"
    ).format(model_name)

    if with_spheno:
        towrite += dumb_indent(4, (
                "// =========================\n"
                "// {0} spectrum (from SARAH-generated SPheno)\n"
                "//\n"
                "#define FUNCTION get_{0}_spectrum_SPheno\n"
                "START_FUNCTION(Spectrum)\n"
                "ALLOW_MODELS({0})\n"
                "DEPENDENCY(SMINPUTS, SMInputs)\n"
                "BACKEND_REQ({0}_spectrum, (libSPheno{0}), int, (Spectrum&, const Finputs&) )\n"
                "BACKEND_OPTION((SARAHSPheno_{0}, {1}), (libSPheno{0}))\n"
                "#undef FUNCTION\n"
                "\n"
        ).format(model_name, SPHENO_VERSION))
    # If we want to make a simple container spectrum only.
    else:
        towrite += (
            "    // Create simple object from SMInputs & new params.\n"
            "    #define FUNCTION get_{0}_spectrum\n"
            "    START_FUNCTION(Spectrum)\n"
            "    DEPENDENCY(SMINPUTS, SMInputs)\n"
        ).format(model_name)

        # Add a conditional Higgs dependency if necessary
        if add_higgs:
            towrite += (
                    "    ALLOW_MODEL_DEPENDENCE(StandardModel_Higgs, {0})\n"
                    "    MODEL_GROUP(higgs, (StandardModel_Higgs))\n"
                    "    MODEL_GROUP({1}, ({0}))\n"
                    "    ALLOW_MODEL_COMBINATION(higgs, {1})\n"
            ).format(model_name, model_name + "_group")
        else:
            towrite += "    ALLOW_MODELS({0})\n"

        towrite +=(
                "    #undef FUNCTION\n"
                "\n"
        )

    towrite += (
        "    // Map for Spectrum, for printing.\n"
        "    #define FUNCTION get_{0}_spectrum_as_map\n"
        "    START_FUNCTION(map_str_dbl)\n"
        "    DEPENDENCY({0}_spectrum, Spectrum)\n"
        "    ALLOW_MODELS({0})\n"
        "    #undef FUNCTION\n"
        "\n"
        "  #undef CAPABILITY\n"
        "\n"
    ).format(model_name)


    if with_spheno and len(higgses) > 1:
        # Go through all Higgses and add the dependencies to the rollcall 
        # for known higgses.
        # If the user has extra higgses you will need to add additional 
        # deps on these within DecayBit and the particle database manually
        higgscontent = ""

        # h0_2
        if 35 in higgses:
        	higgscontent += ( "  DEPENDENCY(Reference_SM_other_Higgs_decay_rates, "
        					  "DecayTable::Entry)\n"  )

        # h0_3 -- we need h0_2 as well
        if 45 in higgses:
        	if 35 not in higgses:
        		raise GumError(("h0_3 in spectrum but not h0_2. Please change "
        						"your PDG code from 45->35."))
        	else:
        		higgscontent += ( "  DEPENDENCY(Reference_SM_h0_3_decay_rates, "
             					  "DecayTable::Entry)\n"  )

        # A0_1
        if 36 in higgses: 
        	higgscontent += ( "  DEPENDENCY(Reference_SM_A0_decay_rates, "
        					  "DecayTable::Entry)\n"  )

        # A0_2 -- we need A0_1 as well
        if 46 in higgses:
        	if 36 not in higgses:
        		raise GumError(("A0_2 in spectrum but not A0_1. Please change "
        						"your PDG code from 46->36."))
        	else:
        		higgscontent += ( "  DEPENDENCY(Reference_SM_A0_2_decay_rates, "
             					  "DecayTable::Entry)\n"  )

        # HiggsBounds output
        towrite += dumb_indent(2, (
                "// Generalised Higgs couplings\n"
                "#define CAPABILITY Higgs_Couplings\n"
                "\n"
                "  // From partial widths\n"
                "  #define FUNCTION {0}_higgs_couplings_SPheno\n"
                "  START_FUNCTION(HiggsCouplingsTable)\n"
                "  DEPENDENCY({0}_spectrum, Spectrum)\n"
                "  // SM rates.\n"
                "  DEPENDENCY(Reference_SM_Higgs_decay_rates, DecayTable::Entry)\n"
                "{1}"
                "  // DecayTable to provide us with the rest of the decays\n"
                "  DEPENDENCY(decay_rates, DecayTable)\n"
                "  BACKEND_REQ({0}_HiggsCouplingsTable, (libSPheno{0}), int, (const Spectrum&, HiggsCouplingsTable&, const Finputs&) )\n"
                "  ALLOW_MODELS({0})\n"
                "  #undef FUNCTION\n"
                "\n"
                "#undef CAPABILITY\n"
                "\n"
        ).format(model_name, higgscontent))

    towrite += "#endif\n"
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
            "double C = alpha_em * pi / (sminputs.GF * pow(2,0.5));\n"
            "double sinW2 = 0.5 - pow( 0.25 - C/pow(sminputs.mZ,2) , 0.5);\n"
            "double cosW2 = 0.5 + pow( 0.25 - C/pow(sminputs.mZ,2) , 0.5);\n"
            "double e = pow( 4*pi*( alpha_em ),0.5);\n"
            "double vev = 1. / sqrt(sqrt(2.)*sminputs.GF);\n"
            "\n"
            "// Gauge couplings\n"
            "{0}.vev = vev;\n"    
            "{0}.g1 = sqrt(5/3) * e / sqrt(cosW2);\n"    
            "{0}.g2 = e / sqrt(sinW2);\n"    
            "{0}.g3 = pow( 4*pi*( sminputs.alphaS ),0.5);\n" 
            "\n"
            "// Yukawas\n"
            "double sqrt2v = pow(2.0,0.5)/vev;\n"
            "{0}.Yu[0][0] = sqrt2v * sminputs.mU;\n"
            "{0}.Yu[1][1] = sqrt2v * sminputs.mCmC;\n"
            "{0}.Yu[2][2] = sqrt2v * sminputs.mT;\n"
            "{0}.Ye[0][0] = sqrt2v * sminputs.mE;\n"
            "{0}.Ye[1][1] = sqrt2v * sminputs.mMu;\n"
            "{0}.Ye[2][2] = sqrt2v * sminputs.mTau;\n"
            "{0}.Yd[0][0] = sqrt2v * sminputs.mD;\n"
            "{0}.Yd[1][1] = sqrt2v * sminputs.mS;\n"
            "{0}.Yd[2][2] = sqrt2v * sminputs.mBmB;\n"
    ).format(model)
      
    return towrite
    
"""
SPHENO-SPECTRA
"""

def write_spheno_spectrum_src(model_name, is_susy):
    """
    Writes SpecBit/src/SpecBit_<MODEL>.cpp
    for a Spectrum object interface to SPheno.
    """

    towrite = (
            "void get_{0}_spectrum_SPheno(Spectrum& spectrum)\n"
            "{{\n"
            "namespace myPipe = Pipes::get_{0}_spectrum_SPheno;\n"
            "const SMInputs &sminputs = *myPipe::Dep::SMINPUTS;\n"
            "\n"
            "// Set up the input structure\n"
            "Finputs inputs;\n"
            "inputs.sminputs = sminputs;\n"
            "inputs.param = myPipe::Param;\n"
            "inputs.options = myPipe::runOptions;\n"
            "\n"
            "// Retrieve any mass cuts\n"
            "static Spectrum::cuts_info mass_cuts = Spectrum::"
            "retrieve_mass_cuts(myPipe::runOptions);\n"
            "\n"
            "// Get the spectrum from the Backend\n"
            "myPipe::BEreq::SARAHSPheno_{0}_spectrum(spectrum, inputs);\n"
            "\n"
            "// Drop SLHA files if requested\n"
            "spectrum.drop_SLHAs_if_requested(myPipe::runOptions, "
            "\"GAMBIT_unimproved_spectrum\");\n"
    ).format(model_name)

    # If it's a susy model, GAMBIT only supports neutralino LSP (for now)
    if (is_susy):
        towrite += (
            "\n"
            "// Only allow neutralino LSPs.\n"
            "if (not has_neutralino_LSP(spectrum)) invalid_point().raise("
            "\"Neutralino is not LSP.\");\n"
        )

    towrite += "}\n"
    
    return towrite


def write_spheno_higgsbounds_interface(model_name, gambit_pdgs, 
                                       neutral_higgses, charged_higgses):
    """
    Writes the HiggsBounds (via SPheno) interface for 
    SpecBit/src/SpecBit_<NewModel>.cpp
    """

    towrite = (
            "\n"
            "/// Put together the Higgs couplings for the {0}, from SPheno\n"
            "void {0}_higgs_couplings_SPheno(HiggsCouplingsTable &result)\n"
            "{{\n"
            "namespace myPipe = Pipes::{0}_higgs_couplings_SPheno;\n"
            "\n"
            "// Retrieve spectrum contents\n"
            "const Spectrum& spec = *myPipe::Dep::{0}_spectrum;\n"
            "const SubSpectrum& he = spec.get_HE();\n"
            "const SMInputs &sminputs = spec.get_SMInputs();\n"
            "\n"
            "const DecayTable* tbl = &(*myPipe::Dep::decay_rates);\n"
            "\n"
            "// Set up the input structure for SPheno\n"
            "Finputs inputs;\n"
            "inputs.sminputs = sminputs;\n"
            "inputs.param = myPipe::Param;\n"
            "inputs.options = myPipe::runOptions;\n"
            "\n"
            "// Set up neutral Higgses\n"
            "static const std::vector<str> sHneut = initVector<str>("
    ).format(model_name)

    # Get the names of all neutral Higgses
    entry = []
    for higgs in neutral_higgses:
        entry.append("\""+pdg_to_particle(higgs, gambit_pdgs)+"\"")

    # Sort the higgses in numerical order - with the neutral ones first
    entry = sorted(entry, key=str.swapcase)

    towrite += ", ".join(entry) + ");\n"
    towrite += (
            "result.set_n_neutral_higgs({0});\n"
            "\n"
            "// Set the CP of the Higgs states. Note that this would\n"
            "// need to be more sophisticated to deal with complex models.\n"
    ).format(len(entry))

    # Assign the CP of each (neutral) Higgs...
    for i in range(len(entry)):
        if entry[i].strip("\"").startswith('h'):
            towrite += (
                "result.CP[{0}] = 1.;  // {1}\n"
            ).format(i, entry[i])
        elif entry[i].strip("\"").startswith('A'):
            towrite += (
                "result.CP[{0}] = -1.; // {1}\n"
            ).format(i, entry[i])
        else:
            raise GumError(("Neutral Higgs with name {} "
                            "found -- gum expects neutral Higgses "
                            "to be named h0_[n] or A0_[n]."
                            ).format(entry[i]))

    # Same deal with the charged Higgses
    towrite += (
            "\n"
            "// Set up charged Higgses\n"
            "static const std::vector<str> sHchar = initVector<str>("
    )

    entry = []
    for higgs in charged_higgses:
        entry.append("\""+pdg_to_particle(higgs, gambit_pdgs)+"\"")

    # Sort the higgses in numerical order - with the neutral ones first
    entry = sorted(entry, key=str.swapcase)

    towrite += ", ".join(entry) + ");\n"
    towrite += (
            "result.set_n_charged_higgs({0});\n"
            "\n"
    ).format(len(entry))

    # We don't need to add the CP for charged Higgses.

    # Get the most 'SM-like' Higgs, needed for HiggsBounds
    # This uses GAMBIT built in functions, only up to NMSSM-like
    # (as of 4/11/19)
    # If the user wishes to include more Higgses, they will need to amend
    # the GAMBIT function SMlike_higgs_PDG_code found in
    # Elements/src/smlike_higgs.cpp
    towrite += (
            "\n"
            "// Work out which SM values correspond to which Higgs\n"
            "int SMlike_higgs = SMlike_higgs_PDG_code(he);\n"
            "\n"
    )

    # Get number of CP even and odd higgses
    nEvenH = len ([pdg_to_particle(x, gambit_pdgs) for x in neutral_higgses if 
                   pdg_to_particle(x, gambit_pdgs).startswith('h')])
    nOddH = len ([pdg_to_particle(x, gambit_pdgs) for x in neutral_higgses if 
                  pdg_to_particle(x, gambit_pdgs).startswith('A')])

    # MSSM-like
    if nEvenH == 2:
        towrite += (
                "// Work out which SM values correspond to which Higgs\n"
                "int higgs = (SMlike_higgs_PDG_code(spec) == 25 ? 0 : 1);\n"
                "int other_higgs = (higgs == 0 ? 1 : 0);\n"
                "\n"
        )

    # NMSSM-like
    elif nEvenH == 3:
        towrite += (
                "int higgs;\n"
                "int other_higgs;\n"
                "int yet_another_higgs;\n"
                "if      (SMlike_higgs == 25) { higgs = 0; other_higgs = 1; "
                "yet_another_higgs = 2; }\n"
                "else if (SMlike_higgs == 35) { higgs = 1; other_higgs = 0; "
                "yet_another_higgs = 2; }\n"
                "else if (SMlike_higgs == 45) { higgs = 2; other_higgs = 0; "
                "yet_another_higgs = 1; }\n"
                "\n"
        )

    # nEvenH > 3 
    else:
        raise GumError(("GAMBIT (and by induction, GUM) does not currently support "
                        "more than 3 CP even Higgses; your model has {}."
                        ).format(nEvenH))

    # Set the SM decays for all Higgses, finally
    towrite += ( 
            "// Set the standard model decays\n"
            "result.set_neutral_decays_SM(higgs, sHneut[higgs], "
            "*myPipe::Dep::Reference_SM_Higgs_decay_rates);\n"
    )
    # h0_2
    if nEvenH > 1:
        towrite += (
                "result.set_neutral_decays_SM(other_higgs, sHneut[other_higgs], "
                "*myPipe::Dep::Reference_SM_other_Higgs_decay_rates);\n"
        )
    # h0_3
    if nEvenH > 2:
        towrite += (
            "result.set_neutral_decays_SM(yet_another_higgs, sHneut[yet_another_higgs], *myPipe::Dep::Reference_SM_h0_3_decay_rates);\n"
        )
    # A0_1
    if nOddH > 0:
        towrite += (
                "result.set_neutral_decays_SM({0}, sHneut[{0}], "
                "*myPipe::Dep::Reference_SM_A0_decay_rates);\n"
        ).format(nEvenH)
    # A0_2
    if nOddH > 1:
        towrite += (
                "result.set_neutral_decays_SM({0}, sHneut[{0}], "
                "*myPipe::Dep::Reference_SM_A0_2_decay_rates);\n"
        ).format(nEvenH+1)

    # Not yet supported
    if nOddH > 2:
        raise GumError(("GAMBIT (and by induction, GUM) does not currently support "
                        "more than 2 CP odd Higgses; your model has {}."
                        ).format(nOddH))
    towrite += (
            "\n"
            "// Set the Higgs sector decays from the DecayTable\n"
            "result.set_neutral_decays(higgs, sHneut[higgs], tbl->at(\"h0_1\"));\n"
    )
    if nEvenH > 1:
        towrite += (
                "result.set_neutral_decays(other_higgs, sHneut[other_higgs], "
                "tbl->at(\"h0_2\"));\n"
        )
    if nEvenH > 2:
        towrite += (
                "result.set_neutral_decays(yet_another_higgs, "
                "sHneut[yet_another_higgs], tbl->at(\"h0_3\"));\n"
        )
    if nOddH > 0:
        towrite += (
                "result.set_neutral_decays({0}, sHneut[{0}], tbl->at(\"A0_1\"));\n"
        ).format(nEvenH)
    if nOddH > 1:
        towrite += (
                "result.set_neutral_decays({0}, sHneut[{0}], tbl->at(\"A0_2\"));\n"
        ).format(nEvenH+1)

    towrite += "\n//Charged Higgses\n"

    # Charged higgses
    for i in range(len(charged_higgses)):
        towrite += (
            "result.set_charged_decays({0}, \"{1}\", tbl->at(\"{1}\"));\n"
        ).format(i, pdg_to_particle(charged_higgses[i], gambit_pdgs))

    # Add the top and wrap it up
    towrite += (
            "\n"
            "// Add t decays since t can decay to light Higgses\n"
            "result.set_t_decays(tbl->at(\"t\"));\n"
            "\n"
            "// Fill HiggsCouplingsTable object from SPheno backend\n"
            "// This fills the effective couplings (C_XX2)\n"
            "myPipe::BEreq::{0}_HiggsCouplingsTable(spec, result, inputs);\n"
            "\n"
            "// The SPheno frontend provides the invisible width for each Higgs, however this requires\n"
            "// loads of additional function calls. Just use the helper function instead.\n"
            "result.invisibles = get_invisibles(he);\n"
            # TODO add invisible particles somehow??
            "}}\n\n"
    ).format(model_name)

    return towrite