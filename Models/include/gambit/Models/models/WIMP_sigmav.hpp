//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  WIMP_sigmav model declarations. 
///  A simple parameterisation of generic WIMP
///  self-annihilation to two-body SM final states.
///  Basic s + p wave expansion, i.e.
///  sigma v = A + B v^2
///  for each channel
///
///  *********************************************
///
///  Authors (add name and date if you modify):
///   
///  \author Ben Farmer
///          (b.farmer@imperial.ac.uk)
///  \date 2019 Aug
///
///  *********************************************


#ifndef __WIMP_sigmav_hpp__
#define __WIMP_sigmav_hpp__

#define MODEL WIMP_sigmav
	START_MODEL
    DEFINEPARS(A_bb, B_bb) 
    DEFINEPARS(A_WW, B_WW) 
    DEFINEPARS(A_cc, B_cc)
    DEFINEPARS(A_tautau, B_tautau)
    DEFINEPARS(A_ZZ, B_ZZ) 
    DEFINEPARS(A_tt, B_tt)
    DEFINEPARS(A_hh, B_hh)

    /// Generic parameterisation of WIMP self-annihilation cross-section to various SM two-body final states
    #define CAPABILITY generic_WIMP_sigmav
    START_CAPABILITY
      #define FUNCTION generic_WIMP_sigmav_from_parameters
      START_FUNCTION(WIMP_annihilation)
      ALLOW_MODEL(WIMP_sigmav)
      #undef FUNCTION
    #undef CAPABILITY

    // Define the module functions
    namespace Gambit {
      namespace Models {
        namespace MODEL {
          void generic_WIMP_sigmav_from_parameters(WIMP_annihilation& result)
          {
              using namespace Pipes::generic_WIMP_sigmav_from_parameters;
              std::vector<std::string> finalstates {"bb", "WW", "cc", "tautau", "ZZ", "tt", "hh"};
              for(auto channel = finalstates.begin(); channel!=finalstates.end(); ++channel)
              {
                 std::string A("A_");
                 std::string B("B_");
                 result.setA(*channel,*Param[A+*channel]);
                 result.setB(*channel,*Param[B+*channel]);
              }
          } 
        }
      }
    }
#undef MODEL

#endif
