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
        inline int CMSSM_spectrum_generator_interface::get_exit_code() const
        {
            return get_BEptr()->get_exit_code();
        }
        
        inline double CMSSM_spectrum_generator_interface::get_reached_precision() const
        {
            return get_BEptr()->get_reached_precision();
        }
        
        inline const flexiblesusy::Spectrum_generator_settings& CMSSM_spectrum_generator_interface::get_settings() const
        {
            return const_cast<flexiblesusy::Abstract_Spectrum_generator_settings&>(const_cast<const Abstract_CMSSM_spectrum_generator_interface<T>*>(get_BEptr())->get_settings__BOSS()).get_init_wref();
        }
        
        inline void CMSSM_spectrum_generator_interface::set_parameter_output_scale(double s)
        {
            get_BEptr()->set_parameter_output_scale(s);
        }
        
        inline void CMSSM_spectrum_generator_interface::set_settings(const flexiblesusy::Spectrum_generator_settings& arg_1)
        {
            get_BEptr()->set_settings__BOSS(*arg_1.get_BEptr());
        }
        
        inline void CMSSM_spectrum_generator_interface::run(const softsusy::QedQcd& arg_1, const flexiblesusy::CMSSM_input_parameters& arg_2)
        {
            get_BEptr()->run__BOSS(*arg_1.get_BEptr(), *arg_2.get_BEptr());
        }
        
        inline void CMSSM_spectrum_generator_interface::write_running_couplings(const ::std::basic_string<char, std::char_traits<char>, std::allocator<char> >& filename, double arg_1, double arg_2) const
        {
            get_BEptr()->write_running_couplings(filename, arg_1, arg_2);
        }
        
        inline void CMSSM_spectrum_generator_interface::write_spectrum(const ::std::basic_string<char, std::char_traits<char>, std::allocator<char> >& filename) const
        {
            get_BEptr()->write_spectrum(filename);
        }
        
        inline void CMSSM_spectrum_generator_interface::write_spectrum() const
        {
            get_BEptr()->write_spectrum__BOSS();
        }
        
        
        // Wrappers for original constructors: 
        inline flexiblesusy::CMSSM_spectrum_generator_interface::CMSSM_spectrum_generator_interface() :
            WrapperBase(__factory0())
        {
            get_BEptr()->set_wptr(this);
            get_BEptr()->set_delete_wrapper(false);
        }
        
        // Special pointer-based constructor: 
        inline flexiblesusy::CMSSM_spectrum_generator_interface::CMSSM_spectrum_generator_interface(flexiblesusy::Abstract_CMSSM_spectrum_generator_interface<T>* in) :
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
        inline flexiblesusy::Abstract_CMSSM_spectrum_generator_interface<T>* flexiblesusy::CMSSM_spectrum_generator_interface::get_BEptr() const
        {
            return dynamic_cast<flexiblesusy::Abstract_CMSSM_spectrum_generator_interface<T>*>(BEptr);
        }
    }
    
}


#include "gambit/Backends/backend_undefs.hpp"

#endif /* __wrapper_CMSSM_spectrum_generator_interface_def_FlexibleSUSY_CMSSM_2_0_1_h__ */
