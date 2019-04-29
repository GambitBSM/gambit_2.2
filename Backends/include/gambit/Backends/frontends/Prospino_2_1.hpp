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

BE_ALLOW_MODELS(MSSM63atQ_mA,MSSM63atMGUT_mA)


// Import functions

// PROSPINO_GB(inlo,isq_ng_in,icoll_in,energy_in,i_error_in,final_state_in,ipart1_in,ipart2_in,isquark1_in,isquark2_in)
// __xx_prospino_subroutine_MOD_prospino_gb
BE_FUNCTION(prospino_gb, void, (Finteger&, Finteger&, Finteger&, Fdouble&, Finteger&, Fstring<2>&, Finteger&, Finteger&, Finteger&, Finteger&, Farray<Fdouble,1,20>&, Farray<Fdouble,0,99>&, Farray<Fdouble,1,2,1,2>&, Farray<Fdouble,1,2,1,2>&, Farray<Fdouble,1,4,1,4>&, Farray<Fdouble,1,2,1,2>&, Farray<Fdouble,1,2,1,2>&, Farray<Fdouble,1,2,1,2>&), "__xx_prospino_subroutine_MOD_prospino_gb", "prospino_gb")

// Convenience functions (registration)
BE_CONV_FUNCTION(run_prospino, std::vector<double>, (const SLHAstruct&, const params_map_type&), "prospino_LHC_xsec")


// Undefine macros to avoid conflict with other backends
#include "gambit/Backends/backend_undefs.hpp"

