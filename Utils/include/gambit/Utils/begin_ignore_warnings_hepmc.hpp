//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  Pragma directives to suppress compiler warnings
///  coming from including HepMC library headers.
///
///  *********************************************
///
///  Authors (add name and date if you modify):
///
///  \author Anders Kvellestad
///          (anders.kvellestad@fys.uio.no)
///  \date 2021 Jul
///
///  *********************************************


#include "gambit/cmake/cmake_variables.hpp"

#ifdef SUPPRESS_LIBRARY_WARNINGS

  // GCC:
  #ifdef __GNUC__
    // Save diagnostic state
    #pragma GCC diagnostic push 
    // Turn off some warnings
    #pragma GCC diagnostic ignored "-Wunused-parameter"
  #endif

  // Clang:
  #ifdef __clang__
    #ifndef __ICC  // icpc apparently also defines __clang__ so need this check too
      // Save diagnostic state
      #pragma clang diagnostic push 
      // Turn off some warnings
      #pragma clang diagnostic ignored "-Wunused-parameter"
    #endif
  #endif

#endif
