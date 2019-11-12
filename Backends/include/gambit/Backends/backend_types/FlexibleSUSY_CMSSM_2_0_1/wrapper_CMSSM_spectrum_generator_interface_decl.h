#ifndef __wrapper_CMSSM_spectrum_generator_interface_decl_FlexibleSUSY_CMSSM_2_0_1_h__
#define __wrapper_CMSSM_spectrum_generator_interface_decl_FlexibleSUSY_CMSSM_2_0_1_h__

#include <cstddef>
#include "forward_decls_wrapper_classes.h"
#include "gambit/Backends/wrapperbase.hpp"
#include "abstract_CMSSM_spectrum_generator_interface.h"
#include "wrapper_Spectrum_generator_settings_decl.h"
#include "wrapper_QedQcd_decl.h"
#include "wrapper_CMSSM_input_parameters_decl.h"
#include <string>

#include "identification.hpp"

namespace CAT_3(BACKENDNAME,_,SAFE_VERSION)
{
    
    namespace flexiblesusy
    {
        
        template<class T>
        class CMSSM_spectrum_generator_interface : public WrapperBase
        {
                // Member variables: 
            public:
                // -- Static factory pointers: 
        
                // -- Other member variables: 
        
                // Member functions: 
        
                // Wrappers for original constructors: 
            public:
        
                // Special pointer-based constructor: 
                CMSSM_spectrum_generator_interface(flexiblesusy::Abstract_CMSSM_spectrum_generator_interface<T>* in);
        
                // Copy constructor: 
                CMSSM_spectrum_generator_interface(const CMSSM_spectrum_generator_interface& in);
        
                // Assignment operator: 
                CMSSM_spectrum_generator_interface& operator=(const CMSSM_spectrum_generator_interface& in);
        
                // Destructor: 
                ~CMSSM_spectrum_generator_interface();
        
                // Returns correctly casted pointer to Abstract class: 
                flexiblesusy::Abstract_CMSSM_spectrum_generator_interface<T>* get_BEptr() const;
        
        };
    }
    
}


#include "gambit/Backends/backend_undefs.hpp"

#endif /* __wrapper_CMSSM_spectrum_generator_interface_decl_FlexibleSUSY_CMSSM_2_0_1_h__ */
