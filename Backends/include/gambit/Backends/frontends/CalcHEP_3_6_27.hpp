//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  Frontend for CalcHEP Backend
///
///  *********************************************
///
///  Authors (add name and date if you modify):
///
///  \author Sanjay Bloor
///          (sanjay.bloor12@imperial.ac.uk)
///  \date 2017 May, Oct
///        2018 Sep
///
///  \author Tomas Gonzalo
///          (tomas.gonzalo@monash.edu)
///  \date 2020 Dec
///
///  *****************************************

#define BACKENDNAME CalcHEP
#define BACKENDLANG CC
#define VERSION 3.6.27
#define SAFE_VERSION 3_6_27

LOAD_LIBRARY

BE_ALLOW_MODELS(ScalarSingletDM_Z2)

BE_FUNCTION(setModel, int, (char*, int), "setModel", "CH_setModel")
BE_FUNCTION(calcMainFunc, int, (), "calcMainFunc", "CH_calcMainFunc")
BE_FUNCTION(getMEcode, numout*, (int ,int, char*, char*, char*, char*), "getMEcode", "CH_getMEcode")
BE_FUNCTION(passParameters, int, (numout*), "passParameters", "CH_passParameters")
BE_FUNCTION(assignVal, int, (char*, double), "assignVal", "CH_assignVal")
BE_FUNCTION(pMass, double, (char*), "pMass", "CH_pMass")
BE_FUNCTION(pdg2mass, char*, (int), "pdg2mass", "CH_pdg2mass")
BE_FUNCTION(pdg2width, char*, (int), "pdg2width", "CH_pdg2width")

BE_VARIABLE(varNames, char*, "varNames", "CH_varNames")
BE_VARIABLE(VZdecay, int, "VZdecay", "CH_VZdecay")
BE_VARIABLE(VWdecay, int, "VWdecay", "CH_VWdecay")

BE_CONV_FUNCTION(generate_decay_code, numout*, (str, str, std::vector<str>), "CH_generate_decay_code")
BE_CONV_FUNCTION(generate_xsec_code, numout*, (str, std::vector<str>, std::vector<str>), "CH_generate_xsec_code")
BE_CONV_FUNCTION(CH_Decay_Width, double, (str&, str&, std::vector<str>&), "CH_CH_Decay_Width")
BE_CONV_FUNCTION(CH_Sigma_V, double, (str&, std::vector<str>&, std::vector<str>&, double&, const DecayTable&), "CH_Sigma_V")
BE_CONV_FUNCTION(Assign_All_Values, void, (const Spectrum&, std::vector<SpectrumParameter>), "CH_Assign_All_Values")
BE_CONV_FUNCTION(Assign_Widths, void, (const DecayTable&), "CH_Assign_Widths")
BE_CONV_FUNCTION(Assign_Value, void, (char*, double), "CH_Assign_Value")

BE_INI_CONDITIONAL_DEPENDENCY(ScalarSingletDM_Z2_spectrum, Spectrum, ScalarSingletDM_Z2)

// Undefine macros to avoid conflict with other backends
#include "gambit/Backends/backend_undefs.hpp"
