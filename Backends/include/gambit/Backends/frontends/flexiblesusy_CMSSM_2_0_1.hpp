//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  Backend macros for FlexibleSUSY_CMSSM 2.4.0
///
///  *********************************************
///
///  Authors (add name and date if you modify):
///
///  \author Peter Athron
///          (peter.athron@monash.edu)
///  \date 2019 Oct
///
///  \author Tomas Gonzalo
///          (tomas.gonzalo@monash.edu)
///  \date 2019 Oct
///
///  *********************************************

#define BACKENDNAME FlexibleSUSY_CMSSM
#define BACKENDLANG CXX
#define VERSION 2.4.0
#define SAFE_VERSION 2_4_0

// Begin
LOAD_LIBRARY

BE_CONV_FUNCTION(run_FS_Spectrum, void, (Spectrum&, const SpectrumInputs&), "FS_CMSSM_Spectrum")

#include "gambit/Backends/backend_undefs.hpp"
