//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  NREO model declarations. 
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

#define MODEL NREO_scalarDM
	START_MODEL
	DEFINEPARS(m)
	DEFINEPARS(c0_1, c0_2, c0_3, c0_4, c0_5, c0_6, c0_7, c0_8, c0_9, c0_10, c0_11, c0_12, c0_13, c0_14, c0_15) 
	DEFINEPARS(c1_1, c1_2, c1_3, c1_4, c1_5, c1_6, c1_7, c1_8, c1_9, c1_10, c1_11, c1_12, c1_13, c1_14, c1_15)

    #define CAPABILITY WIMP_properties
    START_CAPABILITY
       #define FUNCTION NREO_scalarDM_WIMP_properties
       START_FUNCTION(WIMPprops)
       ALLOW_MODELS(NREO_scalarDM)
       #undef FUNCTION
    #undef CAPABILITY

    // Define the module functions
    namespace Gambit {
      namespace Models {
        namespace MODEL {
          void NREO_scalarDM_WIMP_properties(WIMPprops& result)
          {
              using namespace Pipes::NREO_scalarDM_WIMP_properties;
              result.mass   = *Param["m"];
              result.spinx2 = 0;
              result.sc     = true;
              result.name   = "phi";
          } 
        }
      }
    }
#undef MODEL

#define MODEL NREO_MajoranaDM
	START_MODEL
	DEFINEPARS(m)
	DEFINEPARS(c0_1, c0_2, c0_3, c0_4, c0_5, c0_6, c0_7, c0_8, c0_9, c0_10, c0_11, c0_12, c0_13, c0_14, c0_15) 
	DEFINEPARS(c1_1, c1_2, c1_3, c1_4, c1_5, c1_6, c1_7, c1_8, c1_9, c1_10, c1_11, c1_12, c1_13, c1_14, c1_15)

    #define CAPABILITY WIMP_properties
    START_CAPABILITY
       #define FUNCTION NREO_MajoranaDM_WIMP_properties
       START_FUNCTION(WIMPprops)
       ALLOW_MODELS(NREO_MajoranaDM)
       #undef FUNCTION
    #undef CAPABILITY

    // Define the module functions
    namespace Gambit {
      namespace Models {
        namespace MODEL {
          void NREO_MajoranaDM_WIMP_properties(WIMPprops& result)
          {
              using namespace Pipes::NREO_MajoranaDM_WIMP_properties;
              result.mass   = *Param["m"];
              result.spinx2 = 1;
              result.sc     = true;
              result.name   = "chi";
          } 
        }
      }
    }
#undef MODEL

#define MODEL NREO_DiracDM
	START_MODEL
	DEFINEPARS(m)
	DEFINEPARS(c0_1, c0_2, c0_3, c0_4, c0_5, c0_6, c0_7, c0_8, c0_9, c0_10, c0_11, c0_12, c0_13, c0_14, c0_15) 
	DEFINEPARS(c1_1, c1_2, c1_3, c1_4, c1_5, c1_6, c1_7, c1_8, c1_9, c1_10, c1_11, c1_12, c1_13, c1_14, c1_15)

    #define CAPABILITY WIMP_properties
    START_CAPABILITY
       #define FUNCTION NREO_DiracDM_WIMP_properties
       START_FUNCTION(WIMPprops)
       ALLOW_MODELS(NREO_DiracDM)
       #undef FUNCTION
    #undef CAPABILITY

    // Define the module functions
    namespace Gambit {
      namespace Models {
        namespace MODEL {
          void NREO_DiracDM_WIMP_properties(WIMPprops& result)
          {
              using namespace Pipes::NREO_DiracDM_WIMP_properties;
              result.mass   = *Param["m"];
              result.spinx2 = 1;
              result.sc     = false;
              result.name   = "Dchi";
          } 
        }
      }
    }
#undef MODEL



#endif
