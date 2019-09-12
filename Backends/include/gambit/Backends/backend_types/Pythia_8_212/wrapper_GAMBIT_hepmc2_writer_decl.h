#ifndef __wrapper_GAMBIT_hepmc2_writer_decl_Pythia_8_212_h__
#define __wrapper_GAMBIT_hepmc2_writer_decl_Pythia_8_212_h__

#include <cstddef>
#include "forward_decls_wrapper_classes.h"
#include "gambit/Backends/wrapperbase.hpp"
#include "abstract_GAMBIT_hepmc2_writer.h"
#include "wrapper_Pythia_decl.h"

#include "identification.hpp"

namespace CAT_3(BACKENDNAME,_,SAFE_VERSION)
{
    
    namespace Pythia8
    {
        
        class GAMBIT_hepmc2_writer : public WrapperBase
        {
                // Member variables: 
            public:
                // -- Static factory pointers: 
                static Pythia8::Abstract_GAMBIT_hepmc2_writer* (*__factory0)();
        
                // -- Other member variables: 
        
                // Member functions: 
            public:
                void init(const char* filename_in);
        
                void write_event(Pythia8::Pythia* pythia);
        
        
                // Wrappers for original constructors: 
            public:
                GAMBIT_hepmc2_writer();
        
                // Special pointer-based constructor: 
                GAMBIT_hepmc2_writer(Pythia8::Abstract_GAMBIT_hepmc2_writer* in);
        
                // Copy constructor: 
                GAMBIT_hepmc2_writer(const GAMBIT_hepmc2_writer& in);
        
                // Assignment operator: 
                GAMBIT_hepmc2_writer& operator=(const GAMBIT_hepmc2_writer& in);
        
                // Destructor: 
                ~GAMBIT_hepmc2_writer();
        
                // Returns correctly casted pointer to Abstract class: 
                Pythia8::Abstract_GAMBIT_hepmc2_writer* get_BEptr() const;
        
        };
    }
    
}


#include "gambit/Backends/backend_undefs.hpp"

#endif /* __wrapper_GAMBIT_hepmc2_writer_decl_Pythia_8_212_h__ */
