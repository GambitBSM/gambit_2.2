#ifndef __wrapper_Model_decl_FlexibleSUSY_CMSSM_2_0_1_h__
#define __wrapper_Model_decl_FlexibleSUSY_CMSSM_2_0_1_h__

#include <cstddef>
#include "forward_decls_wrapper_classes.h"
#include "gambit/Backends/wrapperbase.hpp"
#include "abstract_Model.h"
#include <string>
#include <ostream>

#include "identification.hpp"

namespace CAT_3(BACKENDNAME,_,SAFE_VERSION)
{
    
    namespace flexiblesusy
    {
        
        class Model : public WrapperBase
        {
                // Member variables: 
            public:
                // -- Static factory pointers: 
                static flexiblesusy::Abstract_Model* (*__factory0)();
        
                // -- Other member variables: 
        
                // Member functions: 
            public:
                void calculate_spectrum();
        
                void clear_problems();
        
                ::std::basic_string<char, std::char_traits<char>, std::allocator<char> > name() const;
        
                void print(::std::basic_ostream<char, std::char_traits<char> >& out) const;
        
                void print() const;
        
                void run_to(double arg_1, double eps);
        
                void run_to(double arg_1);
        
                void set_precision(double arg_1);
        
        
                // Wrappers for original constructors: 
            public:
                Model();
        
                // Special pointer-based constructor: 
                Model(flexiblesusy::Abstract_Model* in);
        
                // Copy constructor: 
                Model(const Model& in);
        
                // Assignment operator: 
                Model& operator=(const Model& in);
        
                // Destructor: 
                ~Model();
        
                // Returns correctly casted pointer to Abstract class: 
                flexiblesusy::Abstract_Model* get_BEptr() const;
        
        };
    }
    
}


#include "gambit/Backends/backend_undefs.hpp"

#endif /* __wrapper_Model_decl_FlexibleSUSY_CMSSM_2_0_1_h__ */
