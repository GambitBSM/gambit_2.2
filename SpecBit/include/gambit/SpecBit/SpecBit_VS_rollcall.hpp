//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  Rollcall declarations for module functions
///  contained in SpecBit_VS.cpp
///
///  *********************************************
///
///  Authors (add name and date if you modify):
///
///  \author James McKay
///          (j.mckay14@imperial.ac.uk)
///  \date Nov 2015
///
///  \author José Eliel Camargo-Molina
///          (elielcamargomolina@gmail.com)
///  \date Jun 2018++
///
///  \author Sanjay Bloor
///          (sanjay.bloor12@imperial.ac.uk)
///  \date Sep 2019
///
///  *********************************************

#ifndef __SpecBit_VS_rollcall_hpp__
#define __SpecBit_VS_rollcall_hpp__

#include "gambit/SpecBit/SpecBit_types.hpp"

  #define CAPABILITY lnL_EW_vacuum
  START_CAPABILITY

    #define FUNCTION check_EW_stability_ScalarSingletDM_Z3
    START_FUNCTION(double)
    DEPENDENCY(ScalarSingletDM_Z3_spectrum, Spectrum)
    ALLOW_MODEL_DEPENDENCE(ScalarSingletDM_Z3)
    MODEL_GROUP(higgs,   (StandardModel_Higgs_running))
    MODEL_GROUP(singlet, (ScalarSingletDM_Z3_running))
    ALLOW_MODEL_COMBINATION(higgs, singlet)
    #undef FUNCTION

  #undef CAPABILITY

  #define CAPABILITY lnL_high_scale_vacuum
  START_CAPABILITY

    #define FUNCTION lnL_highscale_vacuum_decay_single_field
    START_FUNCTION(double)
    DEPENDENCY(high_scale_vacuum_info, dbl_dbl_bool)
    #undef FUNCTION

  #undef CAPABILITY

  #define CAPABILITY high_scale_vacuum_info
  START_CAPABILITY

    #define FUNCTION find_min_lambda_ScalarSingletDM_Z2
    START_FUNCTION(dbl_dbl_bool)
    DEPENDENCY(SMINPUTS, SMInputs)
    DEPENDENCY(ScalarSingletDM_Z2_spectrum, Spectrum)
    ALLOW_MODEL_DEPENDENCE(StandardModel_Higgs_running, ScalarSingletDM_Z2_running)
    MODEL_GROUP(higgs,   (StandardModel_Higgs_running))
    MODEL_GROUP(singlet, (ScalarSingletDM_Z2_running))
    ALLOW_MODEL_COMBINATION(higgs, singlet)
    #undef FUNCTION

    #define FUNCTION find_min_lambda_ScalarSingletDM_Z3
    START_FUNCTION(dbl_dbl_bool)
    DEPENDENCY(SMINPUTS, SMInputs)
    DEPENDENCY(ScalarSingletDM_Z3_spectrum, Spectrum)
    ALLOW_MODEL_DEPENDENCE(StandardModel_Higgs_running, ScalarSingletDM_Z3_running)
    MODEL_GROUP(higgs,   (StandardModel_Higgs_running))
    MODEL_GROUP(singlet, (ScalarSingletDM_Z3_running))
    ALLOW_MODEL_COMBINATION(higgs, singlet)
    #undef FUNCTION

    #define FUNCTION find_min_lambda_MDM
    START_FUNCTION(dbl_dbl_bool)
    DEPENDENCY(SMINPUTS, SMInputs)
    DEPENDENCY(MDM_spectrum, Spectrum)
    ALLOW_MODEL_DEPENDENCE(StandardModel_Higgs_running, MDM)
    MODEL_GROUP(higgs,   (StandardModel_Higgs_running))
    MODEL_GROUP(mdm, (MDM))
    ALLOW_MODEL_COMBINATION(higgs, mdm)
    #undef FUNCTION

  #undef CAPABILITY

  #define CAPABILITY expected_vacuum_lifetime
    START_CAPABILITY
    #define FUNCTION get_expected_vacuum_lifetime
    START_FUNCTION(double)
    DEPENDENCY(high_scale_vacuum_info, dbl_dbl_bool)
    #undef FUNCTION
  #undef CAPABILITY

  #define CAPABILITY check_perturbativity_to_lambda_min
    START_CAPABILITY
    #define FUNCTION check_perturb_min_lambda
    START_FUNCTION(double)
    DEPENDENCY(high_scale_vacuum_info, dbl_dbl_bool)
    #undef FUNCTION
  #undef CAPABILITY

  #define CAPABILITY lambdaB
    START_CAPABILITY
    #define FUNCTION get_lambdaB
    START_FUNCTION(double)
    DEPENDENCY(high_scale_vacuum_info, dbl_dbl_bool)
    #undef FUNCTION
  #undef CAPABILITY

  /**********************/
  /* VEVACIOUS ROUTINES */
  /**********************/

  /// Model dependent -- just tells vevacious the name and locations of the ini files
  /// for each model, since they might not be just <MODELNAME>.vin, etc.
  #define CAPABILITY vevacious_file_location
  START_CAPABILITY
    #define FUNCTION vevacious_file_location_MSSM
    START_FUNCTION(map_str_str)
    #undef FUNCTION
  #undef CAPABILITY

  // Initialize vevacious with a set of YAML runOptions.
  #define CAPABILITY init_vevacious
  START_CAPABILITY
    #define FUNCTION initialize_vevacious
    START_FUNCTION(std::string)
    DEPENDENCY(vevacious_file_location, map_str_str)
    #undef FUNCTION
  #undef CAPABILITY

  // Function to pass spectra to vevacious (via SLHAea). Model dependent.
  #define CAPABILITY pass_spectrum_to_vevacious
  START_CAPABILITY
    #define FUNCTION pass_MSSM_spectrum_to_vevacious
    START_FUNCTION(const vevacious_1_0::VevaciousPlusPlus::VevaciousPlusPlus*)
    DEPENDENCY(unimproved_MSSM_spectrum, Spectrum)
    DEPENDENCY(init_vevacious, std::string)
    ALLOW_MODEL_DEPENDENCE(MSSM, CMSSM)
  #undef FUNCTION
  #undef CAPABILITY

  /// Function for computing the stability of the scalar potential. Model independent. 
  /// Just works with an initialised vevaciousPlusPlus object.
  #define CAPABILITY check_vacuum_stability
  START_CAPABILITY
    #define FUNCTION check_vacuum_stability_vevacious
    START_FUNCTION(SpecBit::VevaciousResultContainer)
    DEPENDENCY(pass_spectrum_to_vevacious, const vevacious_1_0::VevaciousPlusPlus::VevaciousPlusPlus*)
    DEPENDENCY(init_vevacious, std::string)
  #undef FUNCTION
  #undef CAPABILITY

  // Tunnelling likelihood
  #define CAPABILITY VS_likelihood
  START_CAPABILITY
    #define FUNCTION get_likelihood_VS
      START_FUNCTION(double)
      DEPENDENCY(check_vacuum_stability, SpecBit::VevaciousResultContainer)
    #undef FUNCTION
  #undef CAPABILITY

  // Thermal likelihood
  #define CAPABILITY VS_likelihood_thermal
  START_CAPABILITY
    #define FUNCTION get_likelihood_VS_thermal
      START_FUNCTION(double)
      DEPENDENCY(check_vacuum_stability, SpecBit::VevaciousResultContainer)
    #undef FUNCTION
  #undef CAPABILITY

  #define CAPABILITY VS_StraightPathGoodEnough
  START_CAPABILITY
    #define FUNCTION print_VS_StraightPathGoodEnough
      START_FUNCTION(int)
      DEPENDENCY(check_vacuum_stability, SpecBit::VevaciousResultContainer)
    #undef FUNCTION
  #undef CAPABILITY
  
  #define CAPABILITY VS_StraightPathGoodEnough_Thermal
  START_CAPABILITY
    #define FUNCTION print_VS_StraightPathGoodEnough_Thermal
      START_FUNCTION(int)
      DEPENDENCY(check_vacuum_stability, SpecBit::VevaciousResultContainer)
    #undef FUNCTION
  #undef CAPABILITY

#define CAPABILITY VS_ThresholdAndBounceActions
  START_CAPABILITY
    #define FUNCTION print_VS_ThresholdAndBounceActions
      START_FUNCTION(std::vector<double>)
      DEPENDENCY(check_vacuum_stability, SpecBit::VevaciousResultContainer)
    #undef FUNCTION
  #undef CAPABILITY
  
  #define CAPABILITY VS_ThresholdAndBounceActions_Thermal
  START_CAPABILITY
    #define FUNCTION print_VS_ThresholdAndBounceActions_Thermal
      START_FUNCTION(std::vector<double>)
      DEPENDENCY(check_vacuum_stability, SpecBit::VevaciousResultContainer)
    #undef FUNCTION
  #undef CAPABILITY


#endif

