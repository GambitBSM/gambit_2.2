//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  Utility functions for bibtex files
///
///  *********************************************
///
///  Authors (add name and date if you modify):
///
///  \author Tomas Gonzalo
///          (gonzalo@physik.rwth-aachen.de)
///  \date 2021 Sep
///
///  *********************************************


#ifndef __bibtex_functions_hpp__
#define __bibtex_functions_hpp__

#include <vector>
#include <fstream>

#include "gambit/Utils/util_types.hpp"

namespace Gambit
{

  namespace Utils
  {

    // Get the list of bibtex entries on a bibtex file
    std::vector<str> getBibTeXEntries(str);

  }

}

#endif //defined __bibtex_functions_hpp__
