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
#include "gambit/ColliderBit/hepmc2heputils.hpp"
#include "gambit/Utils/util_functions.hpp"
#include "HepMC/IO_GenEvent.h"

namespace Gambit
{

  namespace ColliderBit
  {

    /// A nested function that reads in HepMC event files and converts them to HEPUtils::Event format
    void getLHEvent(HEPUtils::Event& result)
    {
      using namespace Pipes::getLHEvent;

      result.clear();

      // Get the filename and initialise the LHEF reader
      const static str hepmc_filename = runOptions->getValue<str>("hepmc_filename");
      static bool first = true;
      if (first)
      {
        if (not Utils::file_exists(hepmc_filename)) throw std::runtime_error("LHE file "+hepmc_filename+" not found.  Quitting...");
        first = false;
      }
      static HepMC::IO_GenEvent hepmcio(hepmc_filename, std::ios::in);

      // Don't do anything during special iterations
      if (*Loop::iteration < 0) return;

      // Attempt to read the next HepMC event as a HEPUtils event. If there are no more events, wrap up the loop and skip the rest of this iteration.
      bool event_retrieved = true;
      #pragma omp critical (reading_HepMCEvent)
      {
        std::unique_ptr<const HepMC::GenEvent> ge = in.read_next_event();
        if (ge) get_HEPUtils_event(ge, result);
        else event_retrieved = false;
      }
      if (not event_retrieved) Loop::halt();
    }

  }

}
