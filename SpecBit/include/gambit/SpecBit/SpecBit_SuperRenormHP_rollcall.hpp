//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  Rollcall declarations for module functions
///  contained in SpecBit_SingletDM.cpp
///
///  *********************************************
///
///  Authors (add name and date if you modify):
///
///  \author Inigo Saez Casares
///          (inigo.saez_casares@ens-paris-saclay.fr)
///    \date 2020 March
///
///  *********************************************

#ifndef __SpecBit_SuperRenormHP_hpp__
#define __SpecBit_SuperRenormHP_hpp__

  // Spectrum object for SuperRenormHP model
  #define CAPABILITY SuperRenormHP_spectrum
  START_CAPABILITY

    // Create Spectrum object from SMInputs structs, SM Higgs parameters,
    // and the SuperRenormHP parameters
    #define FUNCTION get_SuperRenormHP_spectrum
    START_FUNCTION(Spectrum)
    DEPENDENCY(SMINPUTS, SMInputs)
    ALLOW_MODEL_DEPENDENCE(StandardModel_Higgs, SuperRenormHP)
    MODEL_GROUP(higgs,   (StandardModel_Higgs))
    MODEL_GROUP(singlet, (SuperRenormHP))
    ALLOW_MODEL_COMBINATION(higgs, singlet)
    #undef FUNCTION

    /* // Convert spectrum into a standard map so that it can be printed */
    /* #define FUNCTION get_ScalarSingletDM_Z2_spectrum_as_map */
    /* START_FUNCTION(map_str_dbl) // Just a string to double map. Can't have commas in macro input */
    /* DEPENDENCY(ScalarSingletDM_Z2_spectrum, Spectrum) */
    /* #undef FUNCTION */

  #undef CAPABILITY


  /* // Generalised Higgs couplings */
  /* #define CAPABILITY Higgs_Couplings */

  /*   #define FUNCTION ScalarSingletDM_higgs_couplings_pwid */
  /*   START_FUNCTION(HiggsCouplingsTable) */
  /*   DEPENDENCY(Reference_SM_Higgs_decay_rates, DecayTable::Entry) */
  /*   DEPENDENCY(Higgs_decay_rates, DecayTable::Entry) */
  /*   MODEL_CONDITIONAL_DEPENDENCY(ScalarSingletDM_Z2_spectrum, Spectrum, ScalarSingletDM_Z2, ScalarSingletDM_Z2_running) */
  /*   MODEL_CONDITIONAL_DEPENDENCY(ScalarSingletDM_Z3_spectrum, Spectrum, ScalarSingletDM_Z3, ScalarSingletDM_Z3_running) */
  /*   ALLOW_MODELS(ScalarSingletDM_Z2, ScalarSingletDM_Z2_running, ScalarSingletDM_Z3, ScalarSingletDM_Z3_running) */
  /*   #undef FUNCTION */

  /* #undef CAPABILITY */

  /* // Find scale at which spectrum becomes non-perturbative */
  /* #define CAPABILITY scale_of_nonperturbativity */
  /* START_CAPABILITY */
/* #if(FS_MODEL_ScalarSingletDM_Z3_IS_BUILT) */
  /*   #define FUNCTION find_non_perturb_scale_ScalarSingletDM_Z3 */
  /*   START_FUNCTION(double) */
  /*   DEPENDENCY(ScalarSingletDM_Z3_spectrum, Spectrum) */
  /*   ALLOW_MODELS(ScalarSingletDM_Z3_running) */
  /*   #undef FUNCTION */
/* #endif */

/* #if(FS_MODEL_ScalarSingletDM_Z2_IS_BUILT) */
  /*   #define FUNCTION find_non_perturb_scale_ScalarSingletDM_Z2 */
  /*   START_FUNCTION(double) */
  /*   DEPENDENCY(ScalarSingletDM_Z2_spectrum, Spectrum) */
  /*   ALLOW_MODELS(ScalarSingletDM_Z2,ScalarSingletDM_Z2_running) */
  /*   #undef FUNCTION */
/* #endif */

  /* #undef CAPABILITY */


#endif
