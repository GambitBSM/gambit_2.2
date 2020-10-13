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
///  \author Iñigo Saez Casares
///          (inigo.saez_casares@ens-paris-saclay.fr)
///  \date 2019 December
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

  const Higgs_Nucleon_coupling_fN couplings = *Dep::get_Higgs_Nucleon_coupling_fN;
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

#define FRIEND DecayingDM_mixture

void MODEL_NAMESPACE::SuperRenormHP_to_DecayingDM_mixture (const ModelParameters &myparams, ModelParameters &friendparams)
{
  USE_MODEL_PIPE(FRIEND) // get pipe for "interpret as friend" function
  logger()<<"Running interpret_as_friend calculations for SuperRenormHP -> DecayingDM_mixture ..."<<EOM;

  double mS = myparams["mS"];
  double theta = myparams["theta"];

  double me = Gambit::m_electron;
  double alphaEM = Gambit::alpha_EM;
  double C = 50./27.;
  double pi = Gambit::pi;
  double vev = 256;

  double gamma_ph = (theta*theta*alphaEM*alphaEM*mS*mS*mS*C*C)/(256.*pi*pi*pi*vev*vev);
  double gamma_e = ( mS >= 2*me ) ? pow(theta,2)*pow(me,2)*mS/(8*pi*pow(vev,2))*pow(1 - 4*pow(me,2)/pow(mS,2), 3./2.) : 0;
  double gamma_tot = gamma_ph + gamma_e;

  double RD = *Dep::DM_relic_density;
  double H0 = *Dep::H0;
  double Omega0_cdm = *Dep::Omega0_cdm;

  const double Mpc_2_km = 3.0857e19; // Mpc to km

  double H0_s = H0/Mpc_2_km; // H0 in 1/s
  double rhoC = 3*pow(H0_s, 2)*pow(Gambit::m_planck, 2)/(8*pi)/Gambit::hbar/pow(Gambit::s2cm, 3); // critical density un Gev/cm^3

  friendparams.setValue("mass", mS);
  friendparams.setValue("lifetime", 1/gamma_tot/Gambit::hbar);
  friendparams.setValue("fraction", RD/Omega0_cdm/rhoC);
  friendparams.setValue("BR_ph", gamma_ph/gamma_tot);
  friendparams.setValue("BR_el", gamma_e/gamma_tot);

}
#undef FRIEND
#undef MODEL
