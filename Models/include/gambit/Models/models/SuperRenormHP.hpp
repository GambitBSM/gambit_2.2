//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  Super Renormalizable Higgs Portal DM
///
///  *********************************************
///
///  Authors (add name and date if you modify)
///
///  \author Iñigo Saez Casares
///  \date 2019 December
///
///  *********************************************

#ifndef __SuperRenormHP_hpp__
#define __SuperRenormHP_hpp__

#include "gambit/Models/models/ModifiedGravityYukawa.hpp"

#define MODEL SuperRenormHP
  START_MODEL
  // Physical units : mass [GeV], theta [dimensionless]
  DEFINEPARS(mS) // mass of the DM scalar particle
  DEFINEPARS(theta) // mixing angle with the SM Higgs boson

  // Friendship with ModifiedGravityYukawa
  INTERPRET_AS_X_FUNCTION(ModifiedGravityYukawa,SuperRenormHP_to_ModifiedGravityYukawa)
  INTERPRET_AS_X_DEPENDENCY(ModifiedGravityYukawa, Higgs_Nucleon_coupling_fN, Higgs_Nucleon_coupling_fN)

  MAP_TO_CAPABILITY(mS, DM_mass)
#undef MODEL

#endif
