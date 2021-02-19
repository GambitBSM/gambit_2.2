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

#include <map>
#include <set>
#include <string>
#include <exception>
#include <vector>
#include <utility>
#include <iostream>

#include "gambit/Utils/util_macros.hpp"
#include "gambit/Logs/log_tags.hpp"

#include "gambit/Elements/suspicious_points.hpp"

#include "gambit/Printers/printermanager.hpp"
#include "gambit/Printers/baseprinter.hpp"


namespace Gambit
{

    /// Constructor
    Suspicious_point_exception::Suspicious_point_exception() {} 

    /// Raise the new suspicious point exception, Print it with a message and a code.
    void Suspicious_point_exception::raise(const std::string& msg,int mycode, bool debug)
    {
      // get the printer pointer
      Printers::BaseBasePrinter& printer = *(get_global_printer_manager()->printerptr);

      printer.print(mycode, "Suspicious Point Code", Printers::get_main_param_id("Suspicious Point Code"), printer.getRank(), Printers::get_point_id());

      if (debug) std::cout << "Point Suspicious: " << msg << std::endl;
    }

}
