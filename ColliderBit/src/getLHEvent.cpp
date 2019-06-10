//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  Les Houches event file reader module function
///
///  *********************************************
///
///  Authors (add name and date if you modify):
///
///  \author Pat Scott
///          (p.scott@imperial.ac.uk)
///  \date 2019 May
///
///  *********************************************

#include "gambit/ColliderBit/ColliderBit_eventloop.hpp"
#include "gambit/ColliderBit/lhef2heputils.hpp"
#include "gambit/Utils/util_functions.hpp"

namespace Gambit
{

  namespace ColliderBit
  {

    /// A nested function that reads in Les Houches Event files and converts them to HEPUtils::Event format
    void getLHEvent(HEPUtils::Event& result)
    {
      using namespace Pipes::getLHEvent;

      result.clear();

      // Get the filename and initialise the LHEF reader
      const static str lhef_filename = runOptions->getValue<str>("lhef_filename");
      static bool first = true;
      if (first)
      {
        if (not Utils::file_exists(lhef_filename)) throw std::runtime_error("LHE file "+lhef_filename+" not found.  Quitting...");
        first = false;
      }
      static LHEF::Reader lhe(lhef_filename);

      // Don't do anything during special iterations
      if (*Loop::iteration < 0) return;

      // Attempt to read the next LHE event as a HEPUtils event. If there are no more events, wrap up the loop and skip the rest of this iteration.
      bool event_retrieved = true;
      #pragma omp critical (reading_LHEvent)
      {
        if (lhe.readEvent()) get_HEPUtils_event(lhe, result);
        else event_retrieved = false;
      }
      if (not event_retrieved) Loop::halt();
    }

  }

}