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

#ifndef STANDALONE
    #include "gambit/Printers/printermanager.hpp"
    #include "gambit/Printers/baseprinter.hpp"
#endif

namespace Gambit
{

  /// Gambit suspicious point exception class.
  class Suspicious_point_exception
  {

    public:

      /// Constructor
      Suspicious_point_exception() {} 

      /// Raise the new suspicious point exception, Print it with a message and a code.
      void raise(const std::string& msg,int mycode=1, bool debug=false)
      {

#ifndef STANDALONE
        // get the printer pointer
        Printers::BaseBasePrinter& printer = *(get_global_printer_manager()->printerptr);

        int ranksus = printer.getRank();
        printer.print(mycode,   "Suspicious Point Code", Printers::get_main_param_id("Suspicious Point Code"), ranksus, Printers::get_point_id());
#endif

        if (debug) std::cout << "Point Suspicious: " << msg << std::endl;
      }
  };

}

#endif
