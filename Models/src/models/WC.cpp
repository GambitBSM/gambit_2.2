///  GAMBIT: Global and Modular BSM Inference Tool
///  *********************************************
///
///  Model translation functions for the flavour EFT models
///
///  *********************************************
///
///  Authors
///  =======
///
///  (add name and date if you modify)
///
///  \author Pat Scott
///          (pat.scott@uq.edu.au)
///  \date 2022 May
///
///  *********************************************


#include "gambit/Models/model_macros.hpp"
#include "gambit/Models/model_helpers.hpp"
#include "gambit/Logs/logger.hpp"

#include "gambit/Models/models/WC.hpp"

void MODEL_NAMESPACE::WC_to_WC_LUV (const ModelParameters &myP, ModelParameters &targetP)
{
  logger()<<"Running interpret_as_parent calculations for WC --> WC_LUV."<<LogTags::info<<EOM;

  targetP.setValue("Re_DeltaC7_tau", myP("Re_DeltaC7"));
  targetP.setValue("Im_DeltaC7_tau", myP("Im_DeltaC7"));
  targetP.setValue("Re_DeltaC9_tau", myP("Re_DeltaC9"));
  targetP.setValue("Im_DeltaC9_tau", myP("Im_DeltaC9"));
  targetP.setValue("Re_DeltaC10_tau", myP("Re_DeltaC10"));
  targetP.setValue("Im_DeltaC10_tau", myP("Im_DeltaC10"));
  targetP.setValue("Re_DeltaCQ1_tau", myP("Re_DeltaCQ1"));
  targetP.setValue("Im_DeltaCQ1_tau", myP("Im_DeltaCQ1"));
  targetP.setValue("Re_DeltaCQ2_tau", myP("Re_DeltaCQ2"));
  targetP.setValue("Im_DeltaCQ2_tau", myP("Im_DeltaCQ2"));

  targetP.setValue("Re_DeltaC7_mu", myP("Re_DeltaC7"));
  targetP.setValue("Im_DeltaC7_mu", myP("Im_DeltaC7"));
  targetP.setValue("Re_DeltaC9_mu", myP("Re_DeltaC9"));
  targetP.setValue("Im_DeltaC9_mu", myP("Im_DeltaC9"));
  targetP.setValue("Re_DeltaC10_mu", myP("Re_DeltaC10"));
  targetP.setValue("Im_DeltaC10_mu", myP("Im_DeltaC10"));
  targetP.setValue("Re_DeltaCQ1_mu", myP("Re_DeltaCQ1"));
  targetP.setValue("Im_DeltaCQ1_mu", myP("Im_DeltaCQ1"));
  targetP.setValue("Re_DeltaCQ2_mu", myP("Re_DeltaCQ2"));
  targetP.setValue("Im_DeltaCQ2_mu", myP("Im_DeltaCQ2"));

  targetP.setValue("Re_DeltaC7_e", myP("Re_DeltaC7"));
  targetP.setValue("Im_DeltaC7_e", myP("Im_DeltaC7"));
  targetP.setValue("Re_DeltaC9_e", myP("Re_DeltaC9"));
  targetP.setValue("Im_DeltaC9_e", myP("Im_DeltaC9"));
  targetP.setValue("Re_DeltaC10_e", myP("Re_DeltaC10"));
  targetP.setValue("Im_DeltaC10_e", myP("Im_DeltaC10"));
  targetP.setValue("Re_DeltaCQ1_e", myP("Re_DeltaCQ1"));
  targetP.setValue("Im_DeltaCQ1_e", myP("Im_DeltaCQ1"));
  targetP.setValue("Re_DeltaCQ2_e", myP("Re_DeltaCQ2"));
  targetP.setValue("Im_DeltaCQ2_e", myP("Im_DeltaCQ2"));
}

void MODEL_NAMESPACE::WC_to_WC_LR (const ModelParameters& myP, ModelParameters &targetP)
{
  logger()<<"Running interpret_as_friend calculations for WC -> WC_LR..."<<LogTags::info<<EOM;

  // Send all parameter values upstream to matching parameters in friend.
  targetP.setValues(myP);

  targetP.setValue("Re_DeltaC7_Prime", myP("Re_DeltaC7"));
  targetP.setValue("Im_DeltaC7_Prime", myP("Im_DeltaC7"));
  targetP.setValue("Re_DeltaC9_Prime", myP("Re_DeltaC9"));
  targetP.setValue("Im_DeltaC9_Prime", myP("Im_DeltaC9"));
  targetP.setValue("Re_DeltaC10_Prime", myP("Re_DeltaC10"));
  targetP.setValue("Im_DeltaC10_Prime", myP("Im_DeltaC10"));
  targetP.setValue("Re_DeltaCQ1_Prime", myP("Re_DeltaCQ1"));
  targetP.setValue("Im_DeltaCQ1_Prime", myP("Im_DeltaCQ1"));
  targetP.setValue("Re_DeltaCQ2_Prime", myP("Re_DeltaCQ2"));
  targetP.setValue("Im_DeltaCQ2_Prime", myP("Im_DeltaCQ2"));
}

