#ifndef __wrapper_CMSSM_spectrum_generator_def_FlexibleSUSY_CMSSM_2_0_1_h__
#define __wrapper_CMSSM_spectrum_generator_def_FlexibleSUSY_CMSSM_2_0_1_h__

#include <string>
#include "wrapper_QedQcd_decl.h"
#include "wrapper_CMSSM_input_parameters_decl.h"

#include "identification.hpp"

namespace CAT_3(BACKENDNAME,_,SAFE_VERSION)
{
    
    namespace flexiblesusy
    {
        
        // Member functions: 
        // Special pointer-based constructor: 
        inline flexiblesusy::CMSSM_spectrum_generator::CMSSM_spectrum_generator(flexiblesusy::CMSSM_spectrum_generator* in) :
            CMSSM_spectrum_generator_interface(in)
        {
            get_BEptr()->set_wptr(this);
            get_BEptr()->set_delete_wrapper(false);
        }
        
        // Copy constructor: 
        inline flexiblesusy::CMSSM_spectrum_generator::CMSSM_spectrum_generator(const CMSSM_spectrum_generator& in) :
            CMSSM_spectrum_generator_interface(in.get_BEptr()->pointer_copy__BOSS())
        {
            get_BEptr()->set_wptr(this);
            get_BEptr()->set_delete_wrapper(false);
        }
        
        // Assignment operator: 
        inline flexiblesusy::CMSSM_spectrum_generator& CMSSM_spectrum_generator::operator=(const CMSSM_spectrum_generator& in)
        {
            if (this != &in)
            {
                get_BEptr()->pointer_assign__BOSS(in.get_BEptr());
            }
            return *this;
        }
        
        
        // Destructor: 
        inline flexiblesusy::CMSSM_spectrum_generator::~CMSSM_spectrum_generator()
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
        inline flexiblesusy::CMSSM_spectrum_generator* flexiblesusy::CMSSM_spectrum_generator::get_BEptr() const
        {
            return dynamic_cast<flexiblesusy::CMSSM_spectrum_generator*>(BEptr);
        }
    }
    
}


#include "gambit/Backends/backend_undefs.hpp"

#endif /* __wrapper_CMSSM_spectrum_generator_def_FlexibleSUSY_CMSSM_2_0_1_h__ */
