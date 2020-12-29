//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  Suspicious Point exceptions definitions.
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

#include <string>
#include <iostream>
#include <algorithm>
#include <omp.h>

#include "gambit/Utils/mpiwrapper.hpp"
#include "gambit/Utils/util_macros.hpp"
#include "gambit/Utils/exceptions.hpp"
#include "gambit/Utils/standalone_error_handlers.hpp"
#include "gambit/Logs/logger.hpp"

#include "gambit/Elements/suspicious_points.hpp"
#include "gambit/Printers/printermanager.hpp"

namespace Gambit
{
    /// Suspicious point exception class methods.

    /// Constructor
    Suspicious_point_exception::Suspicious_point_exception() : special_exception("GAMBIT suspicious point."), SuspiciousPtID(Printers::get_main_param_id("Suspicious Point Code")) {} 

    /// Raise the new suspicious point exception, i.e. throw it with a message and a code.
    void Suspicious_point_exception::raise(const std::string& msg,int mycode, bool debug)
    {

      // get the printer pointer
      Printers::BaseBasePrinter& printer = *(get_global_printer_manager()->printerptr);

      int ranksus = printer.getRank();
      printer.print(mycode,   "Suspicious Point Code", SuspiciousPtID, ranksus, Printers::get_point_id());
      if (debug) std::cout << "Point Suspicious: " << msg << endl;

    }

}


