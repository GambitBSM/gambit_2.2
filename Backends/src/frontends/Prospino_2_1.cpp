//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  Frontend for Prospino 2.1 backend
///
///  *********************************************
///
///  Authors (add name and date if you modify):
///
///  \author Anders Kvellestad
///  \date 2019 Apr
///
///  *********************************************

#include "gambit/Backends/frontend_macros.hpp"
#include "gambit/Backends/frontends/Prospino_2_1.hpp"

#include "gambit/Elements/spectrum.hpp"
// #include "gambit/Elements/spectrum_factories.hpp"
// #include "gambit/Models/SimpleSpectra/NMSSMSimpleSpec.hpp"

#include "gambit/Utils/version.hpp"

#define BACKEND_DEBUG 0


// Convenience functions (definition)
BE_NAMESPACE
{
  // Convenience function to run SPheno and obtain the spectrum
  std::vector<double> run_prospino(const Spectrum& spectrum)
  {

    std::cout << "DEBUG: run_prospino: Begin..." << std::endl;

    Finteger inlo = 0;
    Finteger isq_ng_in = 1;
    Finteger icoll_in = 2;

    // Call prospino
    prospino_gb(inlo, isq_ng_in, icoll_in);

    std::cout << "DEBUG: run_prospino: ...End" << std::endl;

    // Dummy result
    std::vector<double> xsec_vals;
    xsec_vals.push_back(1.2345);
    xsec_vals.push_back(1.2345 * 0.1);

    return xsec_vals;
  }
}
END_BE_NAMESPACE

// Initialisation function (definition)
BE_INI_FUNCTION{} END_BE_INI_FUNCTION
