#ifndef __abstract_CMSSM_spectrum_generator_interface_FlexibleSUSY_CMSSM_2_0_1_h__
#define __abstract_CMSSM_spectrum_generator_interface_FlexibleSUSY_CMSSM_2_0_1_h__

#include "gambit/Backends/abstractbase.hpp"
#include "forward_decls_abstract_classes.h"
#include "forward_decls_wrapper_classes.h"
#include "wrapper_Spectrum_generator_settings_decl.h"
#include "wrapper_QedQcd_decl.h"
#include "wrapper_CMSSM_input_parameters_decl.h"
#include <string>
#include <cstddef>
#include <iostream>

#include "enum_decl_copies.h"

#include "identification.hpp"

namespace CAT_3(BACKENDNAME,_,SAFE_VERSION)
{
    
    
    namespace flexiblesusy
    {
        template <>
        class Abstract_Two_scale><flexiblesusy::Two_scale> : public virtual AbstractBase
        {
            public:
    
            private:
                CMSSM_spectrum_generator_interface* wptr;
                bool delete_wrapper;
            public:
                CMSSM_spectrum_generator_interface* get_wptr() { return wptr; }
                void set_wptr(CMSSM_spectrum_generator_interface* wptr_in) { wptr = wptr_in; }
                bool get_delete_wrapper() { return delete_wrapper; }
                void set_delete_wrapper(bool del_wrp_in) { delete_wrapper = del_wrp_in; }
    
            public:
                Abstract_Two_scale>()
                {
                    wptr = 0;
                    delete_wrapper = false;
                }
    
                Abstract_Two_scale>(const Abstract_Two_scale>&)
                {
                    wptr = 0;
                    delete_wrapper = false;
                }
    
                Abstract_Two_scale>& operator=(const Abstract_Two_scale>&) { return *this; }
    
                virtual void init_wrapper() =0;
    
                CMSSM_spectrum_generator_interface* get_init_wptr()
                {
                    init_wrapper();
                    return wptr;
                }
    
                CMSSM_spectrum_generator_interface& get_init_wref()
                {
                    init_wrapper();
                    return *wptr;
                }
    
                virtual ~Abstract_Two_scale>() =0;
        };
    }
    
}


#include "gambit/Backends/backend_undefs.hpp"


#endif /* __abstract_CMSSM_spectrum_generator_interface_FlexibleSUSY_CMSSM_2_0_1_h__ */
