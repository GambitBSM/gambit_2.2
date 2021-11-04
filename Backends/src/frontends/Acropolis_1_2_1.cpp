//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  Frontend for Acropolis 1.2.1 backend
///
///  *********************************************
///
///  Authors (add name and date if you modify):
///
///  \author Tomas Gonzalo
///          (gonzalo@physik.rwth-aachen.de)
///  \date 2021 Oct
///
///  *********************************************

#include "gambit/Backends/frontend_macros.hpp"
#include "gambit/Backends/frontends/Acropolis_1_2_1.hpp"

#ifdef HAVE_PYBIND11

  // Convenience functions (definitions)
  BE_NAMESPACE
  {
    using Backends::cast_std_to_np;
    using Backends::cast_np_to_std;

    namespace py = pybind11;

    template<typename T = double>
    using pyArray = typename py::array_t<T>;
    using pyArray_dbl = pyArray<>;

    // Test routine for Acropolis
    double Acropolis_test()
    {

      std::cout << "test" << std::endl;


      return 0.0;
    }
  }
  END_BE_NAMESPACE

#endif

// Initialisation function (definition)
BE_INI_FUNCTION
{
}
END_BE_INI_FUNCTION
