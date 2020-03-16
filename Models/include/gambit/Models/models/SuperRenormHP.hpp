//  GAMBIT: Global and Modular BSM Inference Tool
//  *********************************************
///
///  Super Renormalizable Higgs Portal DM
///
///  *********************************************
///
///  Authors
///  =======
///
///  (add name and date if you modify)
///
///  \author Iñigo Saez Casares
///  \date 2019 December
///
///  
///
///  *********************************************

#ifndef __SuperRenormHP_hpp__
#define __SuperRenormHP_hpp__

#define MODEL SuperRenormHP
  START_MODEL
  DEFINEPARS(mS, theta)

  MAP_TO_CAPABILITY(mS, DM_mass)
#undef MODEL

#endif
