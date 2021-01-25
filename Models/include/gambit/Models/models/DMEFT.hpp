//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
/// Header file for DMEFT
///
///  Authors (add name and date if you modify):    
///       *** Automatically created by GUM ***     
///                                                
///  \author The GAMBIT Collaboration             
///  \date 12:32PM on October 15, 2019
///
///  \author Sanjay Bloor
///         (sanjay.bloor12@imperial.ac.uk)
///  \date Oct 2019
///                                                
///  ********************************************* 

#ifndef __DMEFT_hpp__
#define __DMEFT_hpp__

#define MODEL DMEFT
  START_MODEL

  DEFINEPARS(Lambda, C51, C52, C61, C62, C63, C64, C71, C72)
  DEFINEPARS(C73, C74, C75, C76, C77, C78, C79, C710, mchi)
  DEFINEPARS(mtrunIN)
  #define CAPABILITY WIMP_properties
  START_CAPABILITY
     #define FUNCTION DMEFT_WIMP_properties
     START_FUNCTION(WIMPprops)
     ALLOW_MODELS(DMEFT)
     #undef FUNCTION
  #undef CAPABILITY
  
  // Define the module functions
  namespace Gambit {
    namespace Models {
      namespace MODEL {
        void DMEFT_WIMP_properties(WIMPprops& result)
        {
            using namespace Pipes::DMEFT_WIMP_properties;
            result.mass   = *Param["mchi"];
            result.spinx2 = 1;
            result.sc     = false;
            result.name   = "chi";
        } 
      }
    }
  }

#undef MODEL

#endif
