#include "gambit/Backends/frontend_macros.hpp"
#include "gambit/Backends/frontends/flexiblesusy_CMSSM_2_0_1.hpp"


// Convenience functions (definitions)
BE_NAMESPACE
{
   void run_FS_Spectrum(Spectrum& spec, const SpectrumInputs& Input)
     {
        backend_warning().raise(LOCAL_INFO, "New FS spectrum calculation not implimented yet.")

     }
}
END_BE_NAMESPACE
