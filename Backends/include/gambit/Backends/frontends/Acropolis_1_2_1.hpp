//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  Frontend header for the Acropolis backend
///
///  Compile-time registration of available
///  functions and variables from this backend.
///
///  *********************************************
///
///  Authors (add name and date if you modify):
///
///  \author Tomas Gonzalo
///          (gonzalo@physik.rwth-aachen.de)
///  \date 2021 
///
///  *********************************************


#define BACKENDNAME Acropolis
#define BACKENDLANG Python
#define VERSION 1.2.1
#define SAFE_VERSION 1_2_1
#define REFERENCE Depta:2020mhj

LOAD_LIBRARY

//BE_ALLOW_MODELS(AnnihilatingDM_general, DecayingDM_general)

BE_CONV_FUNCTION(Acropolis_test, double, (), "acropolis_test")


// Undefine macros to avoid conflict with other backends
#include "gambit/Backends/backend_undefs.hpp"
