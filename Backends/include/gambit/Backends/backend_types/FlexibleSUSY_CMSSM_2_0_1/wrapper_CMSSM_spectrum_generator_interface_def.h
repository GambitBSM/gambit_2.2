#ifndef __wrapper_CMSSM_spectrum_generator_interface_def_FlexibleSUSY_CMSSM_2_0_1_h__
#define __wrapper_CMSSM_spectrum_generator_interface_def_FlexibleSUSY_CMSSM_2_0_1_h__

#include "wrapper_Spectrum_generator_settings_decl.h"
#include "wrapper_QedQcd_decl.h"
#include "wrapper_CMSSM_input_parameters_decl.h"
#include <string>

#include "identification.hpp"

namespace CAT_3(BACKENDNAME,_,SAFE_VERSION)
{
    
    namespace flexiblesusy
    {
        
        // Member functions: 
        // Special pointer-based constructor: 
        inline flexiblesusy::CMSSM_spectrum_generator_interface::CMSSM_spectrum_generator_interface(flexiblesusy::CMSSM_spectrum_generator_interface* in) :
            WrapperBase(in)
        {
            get_BEptr()->set_wptr(this);
            get_BEptr()->set_delete_wrapper(false);
        }
        
        // Copy constructor: 
        inline flexiblesusy::CMSSM_spectrum_generator_interface::CMSSM_spectrum_generator_interface(const CMSSM_spectrum_generator_interface& in) :
            WrapperBase(in.get_BEptr()->pointer_copy__BOSS())
        {
            get_BEptr()->set_wptr(this);
            get_BEptr()->set_delete_wrapper(false);
        }
        
        // Assignment operator: 
        inline flexiblesusy::CMSSM_spectrum_generator_interface& CMSSM_spectrum_generator_interface::operator=(const CMSSM_spectrum_generator_interface& in)
        {
            if (this != &in)
            {
                get_BEptr()->pointer_assign__BOSS(in.get_BEptr());
            }
            return *this;
        }
        
        
        // Destructor: 
        inline flexiblesusy::CMSSM_spectrum_generator_interface::~CMSSM_spectrum_generator_interface()
        {
            if (get_BEptr() != 0)
            {
                get_BEptr()->set_delete_wrapper(false);
                if (can_delete_BEptr())
                {
                    delete BEptr;
                    BEptr = 0;
                }
            }
            set_delete_BEptr(false);
        }
        
        // Returns correctly casted pointer to Abstract class: 
        inline flexiblesusy::CMSSM_spectrum_generator_interface* flexiblesusy::CMSSM_spectrum_generator_interface::get_BEptr() const
        {
            return dynamic_cast<flexiblesusy::CMSSM_spectrum_generator_interface*>(BEptr);
        }
    }
    
}


#include "gambit/Backends/backend_undefs.hpp"

#endif /* __wrapper_CMSSM_spectrum_generator_interface_def_FlexibleSUSY_CMSSM_2_0_1_h__ */
