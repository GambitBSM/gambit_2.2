//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  Frontend for Flavio backend
///
///  *********************************************
///
///  Authors (add name and date if you modify):
///
///
/// \author Marcin Chrzaszcz
///         (mchrzasz@cern.ch)
/// \date 2018
///
///  *********************************************

#include <sstream>
#include "gambit/Backends/frontend_macros.hpp"
#include "gambit/Backends/frontends/Flavio.hpp"
//#include "gambit/Backends/backend_types/SuperIso.hpp"

// Convenience functions (definitions)
BE_NAMESPACE
{
  double sm_prediction_CONV(string)
  {
    double result=sm_prediction(string);
    return result;

  }
  
}
END_BE_NAMESPACE  
