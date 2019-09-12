#ifndef __wrapper_GAMBIT_hepmc2_writer_def_Pythia_8_212_h__
#define __wrapper_GAMBIT_hepmc2_writer_def_Pythia_8_212_h__

#include "wrapper_Pythia_decl.h"

#include "identification.hpp"

namespace CAT_3(BACKENDNAME,_,SAFE_VERSION)
{
    
    namespace Pythia8
    {
        
        // Member functions: 
        inline void GAMBIT_hepmc2_writer::init(const char* filename_in)
        {
            get_BEptr()->init(filename_in);
        }
        
        inline void GAMBIT_hepmc2_writer::write_event(Pythia8::Pythia* pythia)
        {
            get_BEptr()->write_event__BOSS((*pythia).get_BEptr());
        }
        
        
        // Wrappers for original constructors: 
        inline Pythia8::GAMBIT_hepmc2_writer::GAMBIT_hepmc2_writer() :
            WrapperBase(__factory0())
        {
            get_BEptr()->set_wptr(this);
            get_BEptr()->set_delete_wrapper(false);
        }
        
        // Special pointer-based constructor: 
        inline Pythia8::GAMBIT_hepmc2_writer::GAMBIT_hepmc2_writer(Pythia8::Abstract_GAMBIT_hepmc2_writer* in) :
            WrapperBase(in)
        {
            get_BEptr()->set_wptr(this);
            get_BEptr()->set_delete_wrapper(false);
        }
        
        // Copy constructor: 
        inline Pythia8::GAMBIT_hepmc2_writer::GAMBIT_hepmc2_writer(const GAMBIT_hepmc2_writer& in) :
            WrapperBase(in.get_BEptr()->pointer_copy__BOSS())
        {
            get_BEptr()->set_wptr(this);
            get_BEptr()->set_delete_wrapper(false);
        }
        
        // Assignment operator: 
        inline Pythia8::GAMBIT_hepmc2_writer& GAMBIT_hepmc2_writer::operator=(const GAMBIT_hepmc2_writer& in)
        {
            if (this != &in)
            {
                get_BEptr()->pointer_assign__BOSS(in.get_BEptr());
            }
            return *this;
        }
        
        
        // Destructor: 
        inline Pythia8::GAMBIT_hepmc2_writer::~GAMBIT_hepmc2_writer()
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
        inline Pythia8::Abstract_GAMBIT_hepmc2_writer* Pythia8::GAMBIT_hepmc2_writer::get_BEptr() const
        {
            return dynamic_cast<Pythia8::Abstract_GAMBIT_hepmc2_writer*>(BEptr);
        }
    }
    
}


#include "gambit/Backends/backend_undefs.hpp"

#endif /* __wrapper_GAMBIT_hepmc2_writer_def_Pythia_8_212_h__ */
