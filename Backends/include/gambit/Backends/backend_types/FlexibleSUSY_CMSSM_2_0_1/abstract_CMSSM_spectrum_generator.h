#ifndef __abstract_CMSSM_spectrum_generator_FlexibleSUSY_CMSSM_2_0_1_h__
#define __abstract_CMSSM_spectrum_generator_FlexibleSUSY_CMSSM_2_0_1_h__

#include "gambit/Backends/abstractbase.hpp"
#include "forward_decls_abstract_classes.h"
#include "forward_decls_wrapper_classes.h"
#include <string>
#include "wrapper_QedQcd_decl.h"
#include "wrapper_CMSSM_input_parameters_decl.h"
#include <cstddef>
#include <iostream>

#include "enum_decl_copies.h"

#include "identification.hpp"

namespace CAT_3(BACKENDNAME,_,SAFE_VERSION)
{
    
    
    namespace flexiblesusy
    {
        template <>
        class Abstract_CMSSM_spectrum_generator : public virtual AbstractBase
        {
            public:
    
                virtual double get_high_scale() const =0;
    
                virtual double get_susy_scale() const =0;
    
                virtual double get_low_scale() const =0;
    
                virtual double get_pole_mass_scale() const =0;
    
                virtual void write_running_couplings(const ::std::basic_string<char, std::char_traits<char>, std::allocator<char> >&) const =0;
    
                virtual void write_running_couplings__BOSS() const =0;
    
            public:
                virtual void pointer_assign__BOSS(Abstract_CMSSM_spectrum_generator*) =0;
                virtual Abstract_CMSSM_spectrum_generator* pointer_copy__BOSS() =0;
    
            private:
                CMSSM_spectrum_generator<>* wptr;
                bool delete_wrapper;
            public:
                CMSSM_spectrum_generator<>* get_wptr() { return wptr; }
                void set_wptr(CMSSM_spectrum_generator<>* wptr_in) { wptr = wptr_in; }
                bool get_delete_wrapper() { return delete_wrapper; }
                void set_delete_wrapper(bool del_wrp_in) { delete_wrapper = del_wrp_in; }
    
            public:
                Abstract_CMSSM_spectrum_generator()
                {
                    wptr = 0;
                    delete_wrapper = false;
                }
    
                Abstract_CMSSM_spectrum_generator(const Abstract_CMSSM_spectrum_generator&)
                {
                    wptr = 0;
                    delete_wrapper = false;
                }
    
                Abstract_CMSSM_spectrum_generator& operator=(const Abstract_CMSSM_spectrum_generator&) { return *this; }
    
                virtual void init_wrapper() =0;
    
                CMSSM_spectrum_generator<>* get_init_wptr()
                {
                    init_wrapper();
                    return wptr;
                }
    
                CMSSM_spectrum_generator<>& get_init_wref()
                {
                    init_wrapper();
                    return *wptr;
                }
    
                virtual ~Abstract_CMSSM_spectrum_generator() =0;
        };
    }
    
}


#include "gambit/Backends/backend_undefs.hpp"


#endif /* __abstract_CMSSM_spectrum_generator_FlexibleSUSY_CMSSM_2_0_1_h__ */
