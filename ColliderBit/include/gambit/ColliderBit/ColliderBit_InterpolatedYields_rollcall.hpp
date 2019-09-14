//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  Rollcall header for ColliderBit interpolated yields
///  functions.
///
///  *********************************************
///
///  Authors (add name and date if you modify):
///
///  \author Martin White
///          (martin.white@adelaide.edu.au)
///  \date 2019 August
///
///  *********************************************

#include "gambit/ColliderBit/ColliderBit_types.hpp"

#pragma once

#define MODULE ColliderBit

  // Yields defined via interpolation of cross-section and efficiency

  /*#define CAPABILITY InterpolatedAnalysisResults
  START_CAPABILITY
    // DMEFT test model yield
    #define FUNCTION DMEFT_results
  START_FUNCTION(AnalysisDataPointers)
    ALLOW_MODELS(DMEFT_WC_test)
    #undef FUNCTION
    #undef CAPABILITY*/

  #define CAPABILITY DMEFT_ColliderLogLikes
  START_CAPABILITY
    #define FUNCTION calc_DMEFT_ColliderLogLike
    START_FUNCTION(double)
    DEPENDENCY(InterpolatedAnalysisResults, AnalysisNumbers)
    BACKEND_REQ_FROM_GROUP(lnlike_marg_poisson, lnlike_marg_poisson_lognormal_error, (), double, (const int&, const double&, const double&, const double&) )
    BACKEND_REQ_FROM_GROUP(lnlike_marg_poisson, lnlike_marg_poisson_gaussian_error, (), double, (const int&, const double&, const double&, const double&) )
    BACKEND_GROUP(lnlike_marg_poisson)
    #undef FUNCTION
  #undef CAPABILITY


#undef MODULE
