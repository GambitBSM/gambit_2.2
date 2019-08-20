//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  NREO model declaration. 
///
///  *********************************************
///
///  Authors (add name and date if you modify):
///   
///  \author Neal Avis Kozar
///  \date 2018 March
///
///  *********************************************


#ifndef __NREO_hpp__
#define __NREO_hpp__

#define MODEL NREO
	START_MODEL
	DEFINEPARS(j, m, sig_v, wimp_sc)
	DEFINEPARS(c0_1, c0_2, c0_3, c0_4, c0_5, c0_6, c0_7, c0_8, c0_9, c0_10, c0_11, c0_12, c0_13, c0_14, c0_15) 
	DEFINEPARS(c1_1, c1_2, c1_3, c1_4, c1_5, c1_6, c1_7, c1_8, c1_9, c1_10, c1_11, c1_12, c1_13, c1_14, c1_15)

	/// To declare this kind of relationship between a parameter 'my_par' and a 
	/// capability 'capability', one adds the following to the declaration of the 
	/// model containing 'my_par': MAP_TO_CAPABILITY(my_par, capability)
	MAP_TO_CAPABILITY(m, mwimp) /// this capability already exists in DarkBit
    // Need "real" module functions for this
	//MAP_TO_CAPABILITY(j, spinwimpx2) /// As does this (2 times WIMP spin, e.g. 1 means spin 1/2)
    //MAP_TO_CAPABILITY(wimp_sc, wimp_sc) /// " " (boolean flag, 1 if WIMP is self-conjugate, 0 otherwise)
	MAP_TO_CAPABILITY(sig_v, sigmav) /// this capability already exists in DarkBit - not sure if I should actually be setting it, testing

    // Declare a couple of module functions, for the module NREO (automatically created by START_MODEL)
    // This works pretty much the same as in a normal module
    QUICK_FUNCTION(NREO, spinwimpx2, NEW_CAPABILITY, NREO_wimp_spin, unsigned int, (NREO))
    QUICK_FUNCTION(NREO, wimp_sc,    NEW_CAPABILITY, NREO_wimp_sc,   unsigned int, (NREO))

    // Define the module functions
    namespace Gambit {
      namespace Models {
        namespace MODEL {

          // 2x WIMP spin
          void NREO_wimp_spin(unsigned int& result)
          {
              using namespace Pipes::NREO_wimp_spin;
              int twoj = std::round(2*(*Param["j"]));
              if(twoj<0)
              {
                  model_error().raise(LOCAL_INFO,"WIMP spin parameter in NREO model is negative! This parameter must be a postive integer or half-integer!");
              }
              // Could also check that the rounding is less than floating point precision, but currently have neglected to do that
              result = (unsigned int)twoj;
          }

          // WIMP self-conjugateness
          void NREO_wimp_sc(unsigned int& result)
          {
              using namespace Pipes::NREO_wimp_sc;
              int sc = std::round(*Param["wimp_sc"]);
              if(sc!=0 and sc!=1)
              {
                  model_error().raise(LOCAL_INFO,"WIMP sc (self-conjugate) parameter in NREO model is invalid! This parameter must be a postive integer or half-integer!integer");
              }
              // Could also check that the rounding is less than floating point precision, but currently have neglected to do that
              result = (unsigned int)sc;
          }

        }
      }
  }


#undef MODEL

#endif
