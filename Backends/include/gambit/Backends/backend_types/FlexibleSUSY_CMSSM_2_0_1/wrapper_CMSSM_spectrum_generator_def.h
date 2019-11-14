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
        inline double CMSSM_spectrum_generator::get_high_scale() const
        {
            return get_BEptr()->get_high_scale();
        }
        
        inline double CMSSM_spectrum_generator::get_susy_scale() const
        {
            return get_BEptr()->get_susy_scale();
        }
        
        inline double CMSSM_spectrum_generator::get_low_scale() const
        {
            return get_BEptr()->get_low_scale();
        }
        
        inline double CMSSM_spectrum_generator::get_pole_mass_scale() const
        {
            return get_BEptr()->get_pole_mass_scale();
        }
        
        inline void CMSSM_spectrum_generator::write_running_couplings(const ::std::basic_string<char, std::char_traits<char>, std::allocator<char> >& filename) const
        {
            get_BEptr()->write_running_couplings(filename);
        }
        
        inline void CMSSM_spectrum_generator::write_running_couplings() const
        {
            get_BEptr()->write_running_couplings__BOSS();
        }
        
        
        // Wrappers for original constructors: 
        inline flexiblesusy::CMSSM_spectrum_generator::CMSSM_spectrum_generator() :
            WrapperBase(__factory0())
        {
            get_BEptr()->set_wptr(this);
            get_BEptr()->set_delete_wrapper(false);
        }
        
        // Special pointer-based constructor: 
        inline flexiblesusy::CMSSM_spectrum_generator::CMSSM_spectrum_generator(flexiblesusy::Abstract_CMSSM_spectrum_generator<>* in) :
            WrapperBase(in)
        {
            get_BEptr()->set_wptr(this);
            get_BEptr()->set_delete_wrapper(false);
        }
        
        // Copy constructor: 
        inline flexiblesusy::CMSSM_spectrum_generator::CMSSM_spectrum_generator(const CMSSM_spectrum_generator& in) :
            WrapperBase(in.get_BEptr()->pointer_copy__BOSS())
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
        inline flexiblesusy::Abstract_CMSSM_spectrum_generator<>* flexiblesusy::CMSSM_spectrum_generator::get_BEptr() const
        {
            return dynamic_cast<flexiblesusy::Abstract_CMSSM_spectrum_generator<>*>(BEptr);
        }
    }
    
}


#include "gambit/Backends/backend_undefs.hpp"

#endif /* __wrapper_CMSSM_spectrum_generator_def_FlexibleSUSY_CMSSM_2_0_1_h__ */
