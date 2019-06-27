//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  Frontend for HepLike backend v1.0
///
///  *********************************************
///
///  Authors (add name and date if you modify):
///
///  \author Jihyun Bhom
///  \date   2019 July
///
///  *********************************************



#define BACKENDNAME HepLike
#define BACKENDLANG CXX
#define VERSION 1.0
#define SAFE_VERSION 1_0

LOAD_LIBRARY
//BE_FUNCTION(Init_param, void, (parameters*), "Init_param", "Init_param")


// Undefine macros to avoid conflict with other backends
#include "gambit/Backends/backend_undefs.hpp"
