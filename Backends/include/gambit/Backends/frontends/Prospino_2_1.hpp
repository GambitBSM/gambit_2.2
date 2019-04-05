//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  Frontend header for the Prospino backend.
///
///  Compile-time registration of available
///  functions and variables from this backend.
///
///  *********************************************
///
///  Authors (add name and date if you modify):
///
///  \author Anders Kvellestad
///          (a.kvellestad@imperial.ac.uk)
///  \date 2019 March
///
///  *********************************************

#define BACKENDNAME Prospino
#define BACKENDLANG Fortran
#define VERSION 2.1
#define SAFE_VERSION 2_1

// Load it
LOAD_LIBRARY


/*
Functions we'll need:

- PROSPINO_CHECK_HIGGS(final_state_in)
- PROSPINO_CHECK_FS(final_state_in,ipart1_in,ipart2_in,lfinal)
- PROSPINO(inlo,isq_ng_in,icoll_in,energy_in,i_error_in,final_state_in,ipart1_in,ipart2_in,isquark1_in,isquark2_in)

and whatever goes on inside 
- PROSPINO_OPEN_CLOSE(0)
- PROSPINO_OPEN_CLOSE(1)

TODO:
- Make a subroutine PROSPINO_GB that does the same as PROSPINO, 
  but without the file input/output. Perhaps turn it into a 
  function to return the result in an array? Or by reference?

- Get rid of all HARD_STOP calls

*/


// Import functions

// PROSPINO_GB(inlo,isq_ng_in,icoll_in,energy_in,i_error_in,final_state_in,ipart1_in,ipart2_in,isquark1_in,isquark2_in)
// __xx_prospino_subroutine_MOD_prospino_gb
BE_FUNCTION(prospino_gb, void, (), "__xx_prospino_subroutine_MOD_prospino_gb", "prospino_gb")

// BE_FUNCTION(nulike_init, void, (const char&, const char&, const char&, const char&, const char&, double&, bool&, bool&), "nulike_init_", "nulike_init")
// BE_FUNCTION(nulike_bounds, void, (const char&, const double&, const double&, nuyield_function_pointer, double&, double&, int&,
//                                   double&, double&, const int&, const double&, const int&, const bool&, const double&, const double&, void*&, const bool&),
//                                   "nulike_bounds", "nubounds")
// BE_FUNCTION(nulike_lnpiln, double, (const int&, const double&, const double&, const double&), "nulike_lnpiln_", "lnlike_marg_poisson_lognormal_error")
// BE_FUNCTION(nulike_lnpin,  double, (const int&, const double&, const double&, const double&), "nulike_lnpin_",  "lnlike_marg_poisson_gaussian_error")

// Arguments for the last two above are:
//  int    nobs   number of observed events
//  double npred1 number of predicted events with no uncertainty
//  double npred2 number of predicted events with an associated prediction uncertainty due to e.g. efficiency error
//  double error  fractional uncertainty on prediction npred2
// Note that the split into npred1 and npred2 is just for distinguishing which part of the
// predicition has the fractional uncertainty associated with it.  If the uncertainty is on
// the entire prediction, set npred1 = 0 and npred2 = total predicted events.

// Initialisation function (definition)
BE_INI_FUNCTION{} END_BE_INI_FUNCTION

// Undefine macros to avoid conflict with other backends
#include "gambit/Backends/backend_undefs.hpp"

