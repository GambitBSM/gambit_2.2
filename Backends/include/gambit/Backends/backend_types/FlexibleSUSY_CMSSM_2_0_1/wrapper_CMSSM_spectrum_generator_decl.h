#ifndef __wrapper_CMSSM_spectrum_generator_decl_FlexibleSUSY_CMSSM_2_0_1_h__
#define __wrapper_CMSSM_spectrum_generator_decl_FlexibleSUSY_CMSSM_2_0_1_h__

#include <cstddef>
#include "forward_decls_wrapper_classes.h"
#include "gambit/Backends/wrapperbase.hpp"
#include "abstract_CMSSM_spectrum_generator.h"
#include "wrapper_CMSSM_spectrum_generator_interface_decl.h"
#include <string>
#include "wrapper_QedQcd_decl.h"
#include "wrapper_CMSSM_input_parameters_decl.h"

#include "identification.hpp"

namespace CAT_3(BACKENDNAME,_,SAFE_VERSION)
{
    
    namespace flexiblesusy
    {
        
        class CMSSM_spectrum_generator : public CMSSM_spectrum_generator_interface
        {
                // Member variables: 
            public:
                // -- Static factory pointers: 
        
                // -- Other member variables: 
        
                // Member functions: 
        
                // Wrappers for original constructors: 
            public:
        
                // Special pointer-based constructor: 
                CMSSM_spectrum_generator(flexiblesusy::CMSSM_spectrum_generator* in);
        
                // Copy constructor: 
                CMSSM_spectrum_generator(const CMSSM_spectrum_generator& in);
        
                // Assignment operator: 
                CMSSM_spectrum_generator& operator=(const CMSSM_spectrum_generator& in);
        
                // Destructor: 
                ~CMSSM_spectrum_generator();
        
                // Returns correctly casted pointer to Abstract class: 
                flexiblesusy::CMSSM_spectrum_generator* get_BEptr() const;
        
        };
    }
    
}


#include "gambit/Backends/backend_undefs.hpp"

#endif /* __wrapper_CMSSM_spectrum_generator_decl_FlexibleSUSY_CMSSM_2_0_1_h__ */
