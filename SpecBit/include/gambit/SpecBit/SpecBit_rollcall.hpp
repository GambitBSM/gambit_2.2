//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  Rollcall header for module SpecBit
///
///  *********************************************
///
///  Authors (add name and date if you modify):
///
///  \author Tomas Gonzalo
///          (tomas.gonzalo@monash.edu)
///    \date 2019 Oct
///
///  *********************************************

//TODO: This is manual now, we'll automate this with a harvester in the future

#ifndef __SpecBit_rollcall_hpp__
#define __SpecBit_rollcall_hpp__

#define MODULE SpecBit
START_MODULE

// Predefine capabilities that appear in more than one model, to avoid ordering problems
  #define CAPABILITY Higgs_Couplings
  START_CAPABILITY
  #undef CAPABILITY

  #define CAPABILITY scale_of_nonperturbativity
  START_CAPABILITY
  #undef CAPABILITY

#include "gambit/SpecBit/models/DiracSingletDM.hpp"
#include "gambit/SpecBit/models/MajoranaSingletDM.hpp"
#include "gambit/SpecBit/models/MDM.hpp"
#include "gambit/SpecBit/models/MSSM.hpp"
#include "gambit/SpecBit/models/ScalarSingletDM.hpp"
#include "gambit/SpecBit/models/SM.hpp"
#include "gambit/SpecBit/models/VectorSingletDM.hpp"

#include "gambit/SpecBit/SpecBit_VS_rollcall.hpp"

#undef MODULE

#endif /* defined(__SpecBit_rollcall_hpp__) */



