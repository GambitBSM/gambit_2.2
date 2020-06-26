//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  Super Renormalizable Higgs Portal DM
///
///  *********************************************
///
///  Authors (add name and date if you modify):
///
///  \author Sebastian Hoof
///          (s.hoof15@imperial.ac.uk)
///  \date 2019 Feb
///
///  *********************************************

#include <cmath>

#include "gambit/Models/model_macros.hpp"
#include "gambit/Models/model_helpers.hpp"
#include "gambit/Logs/logger.hpp"
#include "gambit/Utils/numerical_constants.hpp"
#include "gambit/DarkBit/DarkBit_types.hpp"

#include "gambit/Models/models/SuperRenormHP.hpp"

#define MODEL SuperRenormHP
#define FRIEND ModifiedGravityYukawa

void MODEL_NAMESPACE::SuperRenormHP_to_ModifiedGravityYukawa (const ModelParameters &myparams, ModelParameters &friendparams)
{
  USE_MODEL_PIPE(FRIEND) // get pipe for "interpret as friend" function
  logger()<<"Running interpret_as_friend calculations for SuperRenormHP -> ModifiedGravityYukawa ..."<<EOM;

  const Higgs_Nucleon_coupling_fN couplings = *Dep::Higgs_Nucleon_coupling_fN;
  const double fN = (couplings.proton+couplings.neutron)/2.;

  const double v = 246; // electroweak vev [GeV]
  const double mS = myparams["mS"], theta = myparams["theta"];
  const double f = v/theta/fN;
  const double Mp = Gambit::m_planck; // Planck mass [GeV]
  const double hbar = Gambit::hbar;  // reduced Planck constant [GeV.s]
  const double cs = Gambit::s2cm*1e-2; // speed of light [m/s]

  friendparams.setValue("alpha", pow(Mp/f, 2)/4./Gambit::pi);
  friendparams.setValue("lambda", hbar*cs/mS);
}
#undef FRIEND
#undef MODEL
