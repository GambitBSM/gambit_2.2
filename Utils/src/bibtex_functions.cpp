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
#include <regex>

#include "gambit/Utils/util_types.hpp"
#include "gambit/Utils/bibtex_functions.hpp"
#include "gambit/Utils/stream_overloads.hpp"

namespace Gambit
{

  namespace Utils
  {

    // Get the list of bibtex entries on a bibtex file
    std::vector<str> getBibTeXEntries(str bibtex_file)
    {
      std::vector<str> entries;
      std::ifstream file(bibtex_file);
      str line;

      // Loop over lines
      while(getline(file, line))
      {
        if(line[0] == '@')
        {
          std::smatch entry;
          std::regex_search(line, entry, std::regex("@.*\\{(.*),"));
          entries.push_back(entry[1].str());
        }
      }
      return entries;
    }

  }

}

#endif //defined __bibtex_functions_hpp__
