#ifndef __wrapper_CMSSM_spectrum_generator_decl_FlexibleSUSY_CMSSM_2_0_1_h__
#define __wrapper_CMSSM_spectrum_generator_decl_FlexibleSUSY_CMSSM_2_0_1_h__

#include <cstddef>
#include "forward_decls_wrapper_classes.h"
#include "gambit/Backends/wrapperbase.hpp"
#include "abstract_CMSSM_spectrum_generator.h"
#include <string>
#include "wrapper_QedQcd_decl.h"
#include "wrapper_CMSSM_input_parameters_decl.h"

#include "identification.hpp"

namespace CAT_3(BACKENDNAME,_,SAFE_VERSION)
{
    
    namespace flexiblesusy
    {
        
        template<>
        class CMSSM_spectrum_generator : public WrapperBase
        {
                // Member variables: 
            public:
                // -- Static factory pointers: 
                static flexiblesusy::Abstract_CMSSM_spectrum_generator<>* (*__factory0)();
        
                // -- Other member variables: 
        
                // Member functions: 
            public:
                double get_high_scale() const;
        
                double get_susy_scale() const;
        
                double get_low_scale() const;
        
                double get_pole_mass_scale() const;
        
                void write_running_couplings(const ::std::basic_string<char, std::char_traits<char>, std::allocator<char> >& filename) const;
        
                void write_running_couplings() const;
        
        
                // Wrappers for original constructors: 
            public:
                CMSSM_spectrum_generator();
        
                // Special pointer-based constructor: 
                CMSSM_spectrum_generator(flexiblesusy::Abstract_CMSSM_spectrum_generator<>* in);
        
                // Copy constructor: 
                CMSSM_spectrum_generator(const CMSSM_spectrum_generator& in);
        
                // Assignment operator: 
                CMSSM_spectrum_generator& operator=(const CMSSM_spectrum_generator& in);
        
                // Destructor: 
                ~CMSSM_spectrum_generator();
        
                // Returns correctly casted pointer to Abstract class: 
                flexiblesusy::Abstract_CMSSM_spectrum_generator<>* get_BEptr() const;
        
        };
    }
    
}


#include "gambit/Backends/backend_undefs.hpp"

#endif /* __wrapper_CMSSM_spectrum_generator_decl_FlexibleSUSY_CMSSM_2_0_1_h__ */
