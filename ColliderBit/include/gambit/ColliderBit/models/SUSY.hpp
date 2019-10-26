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

  // Get Monte Carlo event generator
  #define CAPABILITY HardScatteringSim

    #define FUNCTION getPythia
    START_FUNCTION(Py8Collider_defaultversion)
    NEEDS_MANAGER(RunMC, MCLoopInfo)
    NEEDS_CLASSES_FROM(Pythia, default)
    ALLOW_MODELS(MSSM63atQ, MSSM63atMGUT)
    DEPENDENCY(decay_rates, DecayTable)
    DEPENDENCY(MSSM_spectrum, Spectrum)
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
    #undef FUNCTION

  #undef CAPABILITY

  // Get SLHA content from one or more SLHA files
  #define CAPABILITY SLHAFileNameAndContent
  START_CAPABILITY

    // Get the next SLHA filename and content (for model CB_SLHA_file_model)
    #define FUNCTION getNextSLHAFileNameAndContent
    START_FUNCTION(pair_str_SLHAstruct)
    ALLOW_MODELS(CB_SLHA_file_model)
    #undef FUNCTION  
  
    // Read single SLHA file and replace some entries 
    // (for use with models CB_SLHA_simpmod_scan_model and CB_SLHA_scan_model)
    #define FUNCTION getAndReplaceSLHAContent
    START_FUNCTION(pair_str_SLHAstruct)
    ALLOW_MODELS(CB_SLHA_simpmod_scan_model, CB_SLHA_scan_model)
    #undef FUNCTION  

  #undef CAPABILITY


  // Extract SLHA file elements (for model CB_SLHA_file_model)
  #define CAPABILITY SLHAFileElements
  START_CAPABILITY
    #define FUNCTION getSLHAFileElements
    START_FUNCTION(map_str_dbl)
    ALLOW_MODELS(CB_SLHA_file_model, CB_SLHA_simpmod_scan_model, CB_SLHA_scan_model)
    DEPENDENCY(SLHAFileNameAndContent, pair_str_SLHAstruct)
    #undef FUNCTION
  #undef CAPABILITY

#undef MODULE
