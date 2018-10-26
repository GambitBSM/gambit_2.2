//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  MSSM20atQ_lightgravitino translation function definitions.
///
///  Contains the interpret-as-parent translation
///  functions for:
///
///  MSSM20atQ_lightgravitino        --> MSSM25atQ_lightgravitino
///  MSSM20atQ_mA_lightgravitino     --> MSSM25atQ_mA_lightgravitino
///  MSSM20atMGUT_lightgravitino     --> MSSM25atMGUT_lightgravitino
///  MSSM20atMGUT_mA_lightgravitino  --> MSSM25atMGUT_mA_lightgravitino
///  MSSM20atMSUSY_lightgravitino    --> MSSM25atMSUSY_lightgravitino
///  MSSM20atMSUSY_mA_lightgravitino --> MSSM25atMSUSY_mA_lightgravitino
///
///  As well as the interpret-as-friend translation
///  functions for:
///
///  MSSM20atQ_mA_lightgravitino     --> MSSM20atQ_lightgravitino
///  MSSM20atMGUT_mA_lightgravitino  --> MSSM20atMGUT_lightgravitino
///  MSSM20atMSUSY_mA_lightgravitino --> MSSM20atMSUSY_lightgravitino
///
///  and
///
///  MSSM20atMGUT_lightgravitino  --> MSSM20atQ_lightgravitino
///  MSSM20atMSUSY_lightgravitino --> MSSM20atQ_lightgravitino
///
///  *********************************************
///
///  Authors (add name and date if you modify):
///
///  \author Pat Scott
///          (p.scott@imperial.ac.uk)
///  \date 2015 Sep
///  \date 2018 Oct
///
///  \author Ben Farmer
///          (benjamin.farmer@fysik.su.se)
///  \date 2017 Oct
///
///  *********************************************

#include "gambit/Models/model_macros.hpp"
#include "gambit/Models/model_helpers.hpp"
#include "gambit/Logs/logger.hpp"
#include "gambit/Utils/util_functions.hpp"

#include "gambit/Models/models/MSSM63atQ_lightgravitino.hpp" // Contains declaration of MSSM_mA_lightgravitino_to_MSSM_mhud and MSSMatX_lightgravitino_to_MSSMatQ_lightgravitino functions
#include "gambit/Models/models/MSSM20atQ_lightgravitino.hpp"
#include "gambit/Models/models/MSSM20atQ_mA_lightgravitino.hpp"
#include "gambit/Models/models/MSSM20atMGUT_lightgravitino.hpp"
#include "gambit/Models/models/MSSM20atMGUT_mA_lightgravitino.hpp"
#include "gambit/Models/models/MSSM20atMSUSY_lightgravitino.hpp"
#include "gambit/Models/models/MSSM20atMSUSY_mA_lightgravitino.hpp"

#include "gambit/Elements/spectrum.hpp"


// General helper translation function definition
namespace Gambit
{
  void MSSM20atX_to_MSSM25atX(const ModelParameters &myP, ModelParameters &targetP)
  {
     // Send all parameter values upstream to matching parameters in parent.
     // Ignore that some parameters don't exist in the parent, these are set below.
     targetP.setValues(myP,false);

     // RH squark soft masses, gen 1 and 2
     targetP.setValue("mq2_1",  myP["mq2_12"] ); // mq2_11 in MSSM63
     targetP.setValue("mq2_2",  myP["mq2_12"] ); // mq2_22   " "
     // RH slepton soft masses, gen 1 and 2
     targetP.setValue("ml2_1",  myP["ml2_12"] ); // ml2_11 in MSSM63
     targetP.setValue("ml2_2",  myP["ml2_12"] ); // ml2_22   " "
     // LH down-type squark soft masses
     targetP.setValue("md2_1",  myP["md2_12"] ); // ml2_11 in MSSM63
     targetP.setValue("md2_2",  myP["md2_12"] ); // ml2_22   " "
     // LH up-type squark soft masses
     targetP.setValue("mu2_1",  myP["mu2_12"] ); // mu2_11 in MSSM63
     targetP.setValue("mu2_2",  myP["mu2_12"] ); // mu2_22   " "
     // LH charged slepton soft masses
     targetP.setValue("me2_1",  myP["me2_12"] ); // me2_11 in MSSM63
     targetP.setValue("me2_2",  myP["me2_12"] ); // me2_22   " "
     // Done
  }
}

/// @{ Interpret-as-parent function definitions
/// These are particularly repetitive so let's define them with the help of a macro
#define DEFINE_IAPFUNC(PARENT) \
void MODEL_NAMESPACE::CAT_3(MODEL,_to_,PARENT) (const ModelParameters &myP, ModelParameters &targetP) \
{ \
   logger()<<"Running interpret_as_parent calculations for " STRINGIFY(MODEL) " --> " STRINGIFY(PARENT) "..."<<LogTags::info<<EOM; \
   MSSM20atX_to_MSSM25atX(myP, targetP); \
} \

#define MODEL MSSM20atQ_lightgravitino
DEFINE_IAPFUNC(MSSM25atQ_lightgravitino)
#undef MODEL
#define MODEL MSSM20atQ_mA_lightgravitino
DEFINE_IAPFUNC(MSSM25atQ_mA_lightgravitino)
#undef MODEL
#define MODEL MSSM20atMGUT_lightgravitino
DEFINE_IAPFUNC(MSSM25atMGUT_lightgravitino)
#undef MODEL
#define MODEL MSSM20atMGUT_mA_lightgravitino
DEFINE_IAPFUNC(MSSM25atMGUT_mA_lightgravitino)
#undef MODEL
#define MODEL MSSM20atMSUSY_lightgravitino
DEFINE_IAPFUNC(MSSM25atMSUSY_lightgravitino)
#undef MODEL
#define MODEL MSSM20atMSUSY_mA_lightgravitino
DEFINE_IAPFUNC(MSSM25atMSUSY_mA_lightgravitino)
#undef MODEL
/// @}

/// @{ Interpret-as-friend (mA parameterisations to primary parameterisations)
#define MODEL MSSM20atQ_mA_lightgravitino
void MODEL_NAMESPACE::MSSM20atQ_mA_lightgravitino_to_MSSM20atQ_lightgravitino(const ModelParameters &myP, ModelParameters &targetP)
{
   logger()<<"Running interpret_as_X calculations for MSSM20atQ_mA_lightgravitino --> MSSM20atQ_lightgravitino."<<LogTags::info<<EOM;
   USE_MODEL_PIPE(MSSM20atQ_lightgravitino) // Need the pipe for the TARGET model
   const SubSpectrum& HE = Dep::unimproved_MSSM_spectrum->get_HE();
   MSSM_mA_lightgravitino_to_MSSM_mhud(myP, targetP, HE);
}
#undef MODEL

#define MODEL MSSM20atMGUT_mA_lightgravitino
void MODEL_NAMESPACE::MSSM20atMGUT_mA_lightgravitino_to_MSSM20atMGUT_lightgravitino(const ModelParameters &myP, ModelParameters &targetP)
{
   logger()<<"Running interpret_as_X calculations for MSSM20atMGUT_mA_lightgravitino --> MSSM20atMGUT_lightgravitino."<<LogTags::info<<EOM;
   USE_MODEL_PIPE(MSSM20atMGUT_lightgravitino) // Need the pipe for the TARGET model
   const SubSpectrum& HE = Dep::unimproved_MSSM_spectrum->get_HE();
   MSSM_mA_lightgravitino_to_MSSM_mhud(myP, targetP, HE);
}
#undef MODEL

#define MODEL MSSM20atMSUSY_mA_lightgravitino
void MODEL_NAMESPACE::MSSM20atMSUSY_mA_lightgravitino_to_MSSM20atMSUSY_lightgravitino(const ModelParameters &myP, ModelParameters &targetP)
{
   logger()<<"Running interpret_as_X calculations for MSSM20atMSUSY_mA_lightgravitino --> MSSM20atMSUSY_lightgravitino."<<LogTags::info<<EOM;
   USE_MODEL_PIPE(MSSM20atMSUSY_lightgravitino) // Need the pipe for the TARGET model
   const SubSpectrum& HE = Dep::unimproved_MSSM_spectrum->get_HE();
   MSSM_mA_lightgravitino_to_MSSM_mhud(myP, targetP, HE);
}
#undef MODEL
/// @}

/// @{ Interpret-as-friend (MGUT and MSUSY to Q, in primary parameterisation only)
#define MODEL MSSM20atMGUT_lightgravitino
void MODEL_NAMESPACE::MSSM20atMGUT_lightgravitino_to_MSSM20atQ_lightgravitino (const ModelParameters &myP, ModelParameters &targetP)
{
   USE_MODEL_PIPE(MSSM20atQ_lightgravitino) // Need pipe for TARGET model
   logger()<<"Running interpret_as_X calculations for MSSM20atMGUT_lightgravitino --> MSSM20atQ_lightgravitino..."<<LogTags::info<<EOM;
   const SubSpectrum& HE = Dep::unimproved_MSSM_spectrum->get_HE();
   MSSMatX_to_MSSMatQ_lightgravitino(myP, targetP, HE);
}
#undef MODEL

#define MODEL MSSM20atMSUSY_lightgravitino
void MODEL_NAMESPACE::MSSM20atMSUSY_lightgravitino_to_MSSM20atQ_lightgravitino (const ModelParameters &myP, ModelParameters &targetP)
{
   USE_MODEL_PIPE(MSSM20atQ_lightgravitino) // Need pipe for TARGET model
   logger()<<"Running interpret_as_X calculations for MSSM20atMSUSY_lightgravitino --> MSSM20atQ_lightgravitino..."<<LogTags::info<<EOM;
   const SubSpectrum& HE = Dep::unimproved_MSSM_spectrum->get_HE();
   MSSMatX_to_MSSMatQ_lightgravitino(myP, targetP, HE);
}
#undef MODEL
