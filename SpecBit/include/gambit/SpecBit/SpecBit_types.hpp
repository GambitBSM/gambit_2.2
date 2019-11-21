//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  SpecBit specific types.
///
///  *********************************************
///
///  Authors (add name and date if you modify):
///
///  \author Peter Athron
///          (peter.athron@monash.edu)
///  \date 2019 Oct
///
///  \author Tomas Gonzalo
///          (tomas.gonzalo@monash.edu)
///  \date 2019 Oct
///
///  *********************************************

#ifndef __SpecBit_types_hpp__
#define __SpecBit_types_hpp__

#include "gambit/Utils/util_types.hpp"
#include "gambit/Elements/sminputs.hpp"

namespace Gambit
{

  struct SpectrumInputs
  {
    SMInputs sminputs;
    std::map<str, safe_ptr<const double> > param;
    safe_ptr<Options> options;

    SpectrumInputs(SMInputs smi, std::map<str, safe_ptr<const double> > par, safe_ptr<Options> opt):
      sminputs(smi),
      param(par), 
      options(opt)
    {}
  };

}
#endif //defined __SpecBit_types_hpp__
