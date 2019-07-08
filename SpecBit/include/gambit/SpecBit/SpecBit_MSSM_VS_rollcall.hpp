//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  Rollcall declarations for module functions
///  contained in SpecBit_MSSM_VS.cpp
///
///  *********************************************
///  Authors (add name and date if you modify):
///
///  \author José Eliel Camargo-Molina
///           (elielcamargomolina@gmail.com)
///
///  \date Jun 2018
///
///  *********************************************

#ifndef __SpecBit_MSSM_VS_hpp__
#define __SpecBit_MSSM_VS_hpp__

#include "gambit/SpecBit/SpecBit_types.hpp"


  #define CAPABILITY make_vevaciousPlusPlus_inputs
  START_CAPABILITY
    #define FUNCTION make_vevaciousPlusPlus_inputs
    START_FUNCTION(std::string)
    #undef FUNCTION
  #undef CAPABILITY

  #define CAPABILITY check_stability_MSSM
  START_CAPABILITY
    #define FUNCTION check_stability_MSSM
    START_FUNCTION(SpecBit::VevaciousResultContainer)
    DEPENDENCY(unimproved_MSSM_spectrum, Spectrum)
    DEPENDENCY(make_vevaciousPlusPlus_inputs, std::string)
    ALLOW_MODEL_DEPENDENCE(MSSM, CMSSM)
	#undef FUNCTION
  #undef CAPABILITY

  #define CAPABILITY VS_likelihood_MSSM
  START_CAPABILITY
    #define FUNCTION get_likelihood_VS_MSSM
      START_FUNCTION(double)
      DEPENDENCY(check_stability_MSSM, SpecBit::VevaciousResultContainer)
    #undef FUNCTION
  #undef CAPABILITY

  #define CAPABILITY VS_likelihood_MSSM_thermal
  START_CAPABILITY
    #define FUNCTION get_likelihood_VS_MSSM_thermal
      START_FUNCTION(double)
      DEPENDENCY(check_stability_MSSM, SpecBit::VevaciousResultContainer)
    #undef FUNCTION
  #undef CAPABILITY

  #define CAPABILITY VS_StraighPathGoodEnough
  START_CAPABILITY
    #define FUNCTION print_VS_StraighPathGoodEnough
      START_FUNCTION(int)
      DEPENDENCY(check_stability_MSSM, SpecBit::VevaciousResultContainer)
    #undef FUNCTION
  #undef CAPABILITY
  
  #define CAPABILITY VS_StraighPathGoodEnough_Thermal
  START_CAPABILITY
    #define FUNCTION print_VS_StraighPathGoodEnough_Thermal
      START_FUNCTION(int)
      DEPENDENCY(check_stability_MSSM, SpecBit::VevaciousResultContainer)
    #undef FUNCTION
  #undef CAPABILITY

#define CAPABILITY VS_ThresholdAndBounceActions
  START_CAPABILITY
    #define FUNCTION print_VS_ThresholdAndBounceActions
      START_FUNCTION(std::vector<double>)
      DEPENDENCY(check_stability_MSSM, SpecBit::VevaciousResultContainer)
    #undef FUNCTION
  #undef CAPABILITY
  
  #define CAPABILITY VS_ThresholdAndBounceActions_Thermal
  START_CAPABILITY
    #define FUNCTION print_VS_ThresholdAndBounceActions_Thermal
      START_FUNCTION(std::vector<double>)
      DEPENDENCY(check_stability_MSSM, SpecBit::VevaciousResultContainer)
    #undef FUNCTION
  #undef CAPABILITY



#endif

