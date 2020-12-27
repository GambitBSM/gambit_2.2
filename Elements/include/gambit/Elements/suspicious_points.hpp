//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  Suspicious point exception declarations.
///
///  *********************************************
///
///  Authors (add name and date if you modify):
///
///  \author Chris Chang
///          (christopher.chang@uqconnect.edu.au)
///  \date 2020 Nov
///
///  *********************************************

#ifndef __suspicious_points_hpp__
#define __suspicious_points_hpp__

#include <map>
#include <set>
#include <string>
#include <exception>
#include <vector>
#include <utility>

#include "gambit/Utils/util_macros.hpp"
#include "gambit/Logs/log_tags.hpp"

#include "gambit/Printers/baseprinter.hpp"

namespace Gambit
{

  /// Gambit suspicious point exception class.
  class Suspicious_point_exception : public special_exception
  {

    public:

      const int SuspiciousPtID;

      /// Constructor
      Suspicious_point_exception();

      /// Raise the exception, i.e. throw it. The default code is 1.
      virtual void raise(const std::string&, int code=1, bool debug=false);

  };

}

#endif
