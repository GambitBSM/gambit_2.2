///  GAMBIT: Global and Modular BSM Inference Tool
///  *********************************************
///
///  MSSM25 translation function definitions
///
///  Contains translation functions for
///  MSSM25atQ_lightgravitino        --> MSSM30atQ_lightgravitino
///  MSSM25atQ_mA_lightgravitino     --> MSSM30atQ_mA_lightgravitino
///  MSSM25atMGUT_lightgravitino     --> MSSM30atMGUT_lightgravitino
///  MSSM25atMGUT_mA_lightgravitino  --> MSSM30atMGUT_mA_lightgravitino
///  MSSM25atMSUSY_lightgravitino    --> MSSM30atMSUSY_lightgravitino
///  MSSM25atMSUSY_mA_lightgravitino --> MSSM30atMSUSY_mA_lightgravitino
///
///  *********************************************
///
///  Authors
///  =======
///
///  (add name and date if you modify)
///
///  Ben Farmer
///  2017 Oct
///
///  \author Pat Scott
///          (p.scott@imperial.ac.uk)
///  \date 2018 Oct
///
///  *********************************************

#include <string>
#include <vector>

#include "gambit/Models/model_macros.hpp"
#include "gambit/Models/model_helpers.hpp"
#include "gambit/Logs/logger.hpp"
#include "gambit/Utils/util_functions.hpp"

#include "gambit/Models/models/MSSM25atQ_lightgravitino.hpp"
#include "gambit/Models/models/MSSM25atQ_mA_lightgravitino.hpp"
#include "gambit/Models/models/MSSM25atMGUT_lightgravitino.hpp"
#include "gambit/Models/models/MSSM25atMGUT_mA_lightgravitino.hpp"
#include "gambit/Models/models/MSSM25atMSUSY_lightgravitino.hpp"
#include "gambit/Models/models/MSSM25atMSUSY_mA_lightgravitino.hpp"


using namespace Gambit::Utils;

// General helper translation function definition
namespace Gambit {
  void MSSM25atX_to_MSSM30atX(const ModelParameters &myP, ModelParameters &targetP)
  {
     // Copy all the common parameters of MSSM25atQ_lightgravitino into MSSM30atQ_lightgravitino
     targetP.setValues(myP,false);

     // Manually set the parameters which differ
     // slepton trilinear couplings
     // Off-diagonal elements set to zero by parent model
     // First and second generation elements set equal
     targetP.setValue("Ae_1",  myP["Ae_12"] ); // Ae2_11 in MSSM63
     targetP.setValue("Ae_2",  myP["Ae_12"] ); // Ae2_22   " "
     //targetP.setValue("Ae_3",  myP["Ae_3"]  ); // Ae2_33 // Taken care of by common parameter copy

     // down-type trilinear couplings
     // Off-diagonal elements set to zero by parent model
     // First and second generation to zero
     targetP.setValue("Ad_1",  0. );          // Ad2_11 in MSSM63
     targetP.setValue("Ad_2",  0. );          // Ad2_22   " "
     //targetP.setValue("Ad_3",  myP["Ad_3"] ); // Ad2_33 // Taken care of by common parameter copy

     // up-type trilinear couplings
     // Off-diagonal elements set to zero by parent model
     // First and second generation set to zero
     targetP.setValue("Au_1",  0. );          // Au2_11 in MSSM63
     targetP.setValue("Au_2",  0. );          // Au2_22   " "
     // targetP.setValue("Au_3",  myP["Au_3"] ); // Au2_33 // Taken care of by common parameter copy

     // Done
  }
}

/// @{ Interpret-as-parent function definitions
/// These are particularly repetitive so let's define them with the help of a macro
#define DEFINE_IAPFUNC(PARENT) \
void MODEL_NAMESPACE::CAT_3(MODEL,_to_,PARENT) (const ModelParameters &myP, ModelParameters &targetP) \
{ \
   logger()<<"Running interpret_as_parent calculations for " STRINGIFY(MODEL) " --> " STRINGIFY(PARENT) "..."<<LogTags::info<<EOM; \
   MSSM25atX_to_MSSM30atX(myP, targetP); \
} \

#define MODEL MSSM25atQ_lightgravitino
DEFINE_IAPFUNC(MSSM30atQ_lightgravitino)
#undef MODEL
#define MODEL MSSM25atQ_mA_lightgravitino
DEFINE_IAPFUNC(MSSM30atQ_mA_lightgravitino)
#undef MODEL
#define MODEL MSSM25atMGUT_lightgravitino
DEFINE_IAPFUNC(MSSM30atMGUT_lightgravitino)
#undef MODEL
#define MODEL MSSM25atMGUT_mA_lightgravitino
DEFINE_IAPFUNC(MSSM30atMGUT_mA_lightgravitino)
#undef MODEL
#define MODEL MSSM25atMSUSY_lightgravitino
DEFINE_IAPFUNC(MSSM30atMSUSY_lightgravitino)
#undef MODEL
#define MODEL MSSM25atMSUSY_mA_lightgravitino
DEFINE_IAPFUNC(MSSM30atMSUSY_mA_lightgravitino)
#undef MODEL
/// @}

