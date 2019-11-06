//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  Rollcall header for ColliderBit module;
///  SUSY models.
///
///  *********************************************
///
///  Authors (add name and date if you modify):
///
///  \author Abram Krislock
///          (a.m.b.krislock@fys.uio.no)
///
///  \author Aldo Saavedra
///
///  \author Christopher Rogan
///          (christophersrogan@gmail.com)
///  \date 2015 Apr
///
///  \author Pat Scott
///          (p.scott@imperial.ac.uk)
///  \date 2015 Jul
///  \date 2018 Jan
///
///  \author Andy Buckley
///          (andy.buckley@cern.ch)
///  \date 2017 Jun
///
///  *********************************************

#pragma once

#define MODULE ColliderBit

  // SLHAea object with spectrum and decays for a Pythia8 collider
  #define CAPABILITY SpectrumAndDecaysForPythia

    #define FUNCTION getSpectrumAndDecaysForPythia
    START_FUNCTION(SLHAstruct)
    ALLOW_MODELS(MSSM63atQ, MSSM63atMGUT)
    NEEDS_MANAGER(RunMC, MCLoopInfo)      // @todo Why is this needed?
    DEPENDENCY(decay_rates, DecayTable)
    DEPENDENCY(MSSM_spectrum, Spectrum)
    DEPENDENCY(SLHA_pseudonyms, mass_es_pseudonyms)
    #undef FUNCTION

  #undef CAPABILITY


  // Get Monte Carlo event generator
  #define CAPABILITY HardScatteringSim

    #define FUNCTION getPythia
    START_FUNCTION(Py8Collider_defaultversion)
    NEEDS_MANAGER(RunMC, MCLoopInfo)
    NEEDS_CLASSES_FROM(Pythia, default)
    DEPENDENCY(SpectrumAndDecaysForPythia, SLHAstruct)
    #undef FUNCTION

    #define FUNCTION getPythia_SLHA
    START_FUNCTION(Py8Collider_defaultversion)
    NEEDS_MANAGER(RunMC, MCLoopInfo)
    NEEDS_CLASSES_FROM(Pythia, default)
    ALLOW_MODELS(CB_SLHA_file_model, CB_SLHA_simpmod_scan_model, CB_SLHA_scan_model)
    DEPENDENCY(SLHAFileNameAndContent, pair_str_SLHAstruct)
    #undef FUNCTION

    #define FUNCTION getPythiaAsBase
    START_FUNCTION(const BaseCollider*)
    NEEDS_MANAGER(RunMC, MCLoopInfo)
    NEEDS_CLASSES_FROM(Pythia, default)
    DEPENDENCY(HardScatteringSim, Py8Collider_defaultversion)
    #undef FUNCTION

  #undef CAPABILITY


  // Run event generator
  #define CAPABILITY HardScatteringEvent

    #define FUNCTION generateEventPythia
    START_FUNCTION(HEPUtils::Event)
    NEEDS_MANAGER(RunMC, MCLoopInfo)
    NEEDS_CLASSES_FROM(Pythia, default)
    DEPENDENCY(HardScatteringSim, Py8Collider_defaultversion)
    DEPENDENCY(EventWeighterFunction, EventWeighterFunctionType)
    #undef FUNCTION

  #undef CAPABILITY




  /// Cross-sections for weighting events by production process
  /// @{

  /// A map between PID pairs and cross-sections
  #define CAPABILITY PIDPairCrossSectionsMap

    /// Example of provding PIDPairCrossSectionsMap using a Python backend
    /// @todo 1. Replace SLHA1Spectrum dependency with SpectrumAndDecaysForPythia (to ensure same spectrum)
    /// @todo 2. Add a CB utility function that checks if a SLHAstruct is SLHA1 or SLHA2, and use it in this function
    #define FUNCTION getPIDPairCrossSectionsMap_xsecBE_example
    START_FUNCTION(map_PID_pair_PID_pair_xsec)
    NEEDS_MANAGER(RunMC, MCLoopInfo)
    DEPENDENCY(ActivePIDPairs, vec_PID_pair)
    DEPENDENCY(SLHA1Spectrum, SLHAstruct)
    ALLOW_MODELS(MSSM63atQ, MSSM63atMGUT, CB_SLHA_file_model, CB_SLHA_simpmod_scan_model, CB_SLHA_scan_model)
    BACKEND_REQ(xsecBE_example_xsec_fb, (), double, (iipair&, pybind11::dict&, pybind11::dict&))
    BACKEND_REQ(xsecBE_example_xsec_err_fb, (), ddpair, (iipair&, pybind11::dict&, pybind11::dict&))
    BACKEND_REQ(xsecBE_example_set_parameters, (), void, (pybind11::dict&))
    BACKEND_REQ(xsecBE_example_set_flags, (), void, (pybind11::dict&))
    #undef FUNCTION
  
    /// Example of provding PIDPairCrossSectionsMap using a Python backend
    /// @todo 1. Replace SLHA1Spectrum dependency with SpectrumAndDecaysForPythia (to ensure same spectrum)
    /// @todo 2. Add a CB utility function that checks if a SLHAstruct is SLHA1 or SLHA2, and use it in this function
    #define FUNCTION getPIDPairCrossSectionsMap_xsecBE
    START_FUNCTION(map_PID_pair_PID_pair_xsec)
    NEEDS_MANAGER(RunMC, MCLoopInfo)
    DEPENDENCY(ActivePIDPairs, vec_PID_pair)
    DEPENDENCY(SLHA1Spectrum, SLHAstruct)
    ALLOW_MODELS(MSSM63atQ, MSSM63atMGUT, CB_SLHA_file_model, CB_SLHA_simpmod_scan_model, CB_SLHA_scan_model)
    BACKEND_REQ(xsecBE_import_slha_string, (), void, (std::string&))
    BACKEND_REQ(xsecBE_set_parameters, (), void, (pybind11::dict&))
    BACKEND_REQ(xsecBE_get_xsection, (), pybind11::dict, (iipair&))
    #undef FUNCTION

    // #define FUNCTION getProspinoxsec
    // START_FUNCTION(xsec)
    // NEEDS_MANAGER(RunMC, MCLoopInfo)
    // ALLOW_MODELS(MSSM63atQ, MSSM63atMGUT)
    // DEPENDENCY(MSSM_spectrum, Spectrum)
    // BACKEND_REQ(prospino_LHC_xsec, (libprospino), map_str_dbl, (const SLHAstruct&, const param_map_type&, prospino_settings&))
    // #undef FUNCTION

  #undef CAPABILITY
  /// @}


  /// Get SLHA content from one or more SLHA files
  /// @{
  #define CAPABILITY SLHAFileNameAndContent
  START_CAPABILITY

    /// Get the next SLHA filename and content (for model CB_SLHA_file_model)
    #define FUNCTION getNextSLHAFileNameAndContent
    START_FUNCTION(pair_str_SLHAstruct)
    ALLOW_MODELS(CB_SLHA_file_model)
    #undef FUNCTION

    /// Read single SLHA file and replace some entries
    /// (for use with models CB_SLHA_simpmod_scan_model and CB_SLHA_scan_model)
    #define FUNCTION getAndReplaceSLHAContent
    START_FUNCTION(pair_str_SLHAstruct)
    ALLOW_MODELS(CB_SLHA_simpmod_scan_model, CB_SLHA_scan_model)
    #undef FUNCTION

  #undef CAPABILITY
  /// @}


  /// Extract SLHA file elements (for model CB_SLHA_file_model)
  /// @{
  #define CAPABILITY SLHAFileElements
  START_CAPABILITY
    #define FUNCTION getSLHAFileElements
    START_FUNCTION(map_str_dbl)
    ALLOW_MODELS(CB_SLHA_file_model, CB_SLHA_simpmod_scan_model, CB_SLHA_scan_model)
    DEPENDENCY(SLHAFileNameAndContent, pair_str_SLHAstruct)
    #undef FUNCTION
  #undef CAPABILITY
  /// @}


  /// Extract an SLHAstruct with the spectrum
  /// @{
  #define CAPABILITY SLHA1Spectrum
  START_CAPABILITY
    #define FUNCTION getSLHA1Spectrum
    START_FUNCTION(SLHAstruct)
    ALLOW_MODELS(MSSM63atQ, MSSM63atMGUT, CB_SLHA_file_model, CB_SLHA_simpmod_scan_model, CB_SLHA_scan_model)
    MODEL_CONDITIONAL_DEPENDENCY(SLHAFileNameAndContent, pair_str_SLHAstruct, CB_SLHA_file_model, CB_SLHA_simpmod_scan_model, CB_SLHA_scan_model)
    MODEL_CONDITIONAL_DEPENDENCY(MSSM_spectrum, Spectrum, MSSM63atQ, MSSM63atMGUT)
    #undef FUNCTION
  #undef CAPABILITY

  #define CAPABILITY SLHA2Spectrum
  START_CAPABILITY
    #define FUNCTION getSLHA2Spectrum
    START_FUNCTION(SLHAstruct)
    ALLOW_MODELS(MSSM63atQ, MSSM63atMGUT, CB_SLHA_file_model, CB_SLHA_simpmod_scan_model, CB_SLHA_scan_model)
    MODEL_CONDITIONAL_DEPENDENCY(SLHAFileNameAndContent, pair_str_SLHAstruct, CB_SLHA_file_model, CB_SLHA_simpmod_scan_model, CB_SLHA_scan_model)
    MODEL_CONDITIONAL_DEPENDENCY(MSSM_spectrum, Spectrum, MSSM63atQ, MSSM63atMGUT)
    #undef FUNCTION
  #undef CAPABILITY
  /// @}


#undef MODULE
