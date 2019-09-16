//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  HepMC event file reader module function
///
///  *********************************************
///
///  Authors (add name and date if you modify):
///
///  \author Andy Buckley
///          (andy.buckley@cern.ch)
///  \date 2019 June
///
///  \author Pat Scott
///          (p.scott@imperial.ac.uk)
///  \date 2019 June
///
///  \author Anders Kvellestad
///          (a.kvellestad@imperial.ac.uk)
///  \date 2019 June
///
///  \author Tomas Gonzalo
///          (tomas.gonzalo@monash.edu)
///  \date 2019 Sep
///
///  *********************************************

#include "gambit/cmake/cmake_variables.hpp"

#ifndef EXCLUDE_HEPMC

#include "gambit/ColliderBit/ColliderBit_eventloop.hpp"
#include "gambit/ColliderBit/hepmc2heputils.hpp"
#include "gambit/Utils/util_functions.hpp"
#include "HepMC3/ReaderAsciiHepMC2.h"
#include "HepMC3/ReaderAscii.h"

//#define COLLIDERBIT_DEBUG

namespace Gambit
{

  namespace ColliderBit
  {

    /// A nested function that reads in HepMC event files and converts them to HEPUtils::Event format
    void getHepMCEvent(HEPUtils::Event& result)
    {
      using namespace Pipes::getHepMCEvent;

      result.clear();

      // Get the filename and initialise the HepMC reader
      const static str HepMC2_filename = runOptions->getValueOrDef<str>("", "HepMC2_filename");
      static bool HepMC2_ON = (HepMC2_filename != "");
      const static str HepMC3_filename = runOptions->getValueOrDef<str>("", "HepMC3_filename");
      static bool HepMC3_ON = (HepMC3_filename != "");

      static bool first = true;
      if (first)
      {
        if (HepMC2_ON and HepMC3_ON) throw std::runtime_error("Cannot read simultaneously from HepMC2 and HepMC3 files. Quitting...");
        if (HepMC2_ON)
        {
          if (not Utils::file_exists(HepMC2_filename)) throw std::runtime_error("HepMC2 event file "+HepMC2_filename+" not found.  Quitting...");
        }
        else if (HepMC3_ON)
        {
          if (not Utils::file_exists(HepMC3_filename)) throw std::runtime_error("HepMC3 event file "+HepMC3_filename+" not found.  Quitting...");
        }
        else
          throw std::runtime_error("Neither HepMC2 nor HepMC3 event file found. Quitting...");
        first = false;
      }

      if(HepMC2_ON) 
      {
        static HepMC3::ReaderAsciiHepMC2 HepMC2io(HepMC2_filename);

        // Don't do anything during special iterations
        if (*Loop::iteration < 0) return;

        #ifdef COLLIDERBIT_DEBUG
          cout << "Event number: " << *Loop::iteration << endl;
        #endif

        // Attempt to read the next HepMC event as a HEPUtils event. If there are no more events, wrap up the loop and skip the rest of this iteration.
        HepMC3::GenEvent ge(HepMC3::Units::GEV, HepMC3::Units::MM);
        bool event_retrieved = true;
        #pragma omp critical (reading_HepMCEvent)
        {
          event_retrieved = HepMC2io.read_event(ge);
       
          // FIXME This is a temp solution to ensure that the event reading
          //       stops when there are no more events in the HepMC file.
          //       Remove this once bugfix is implemented in HepMC.
          if ((ge.particles().size() == 0) && (ge.vertices().size() == 0)) event_retrieved = false;
        }
        if (not event_retrieved) Loop::halt();

        // Translate to HEPUtils event
        get_HEPUtils_event(ge, result);
 
      }
      else
      {
        static HepMC3::ReaderAscii HepMC3io(HepMC3_filename);

        // Don't do anything during special iterations
        if (*Loop::iteration < 0) return;

        #ifdef COLLIDERBIT_DEBUG
          cout << "Event number: " << *Loop::iteration << endl;
        #endif

        // Attempt to read the next HepMC event as a HEPUtils event. If there are no more events, wrap up the loop and skip the rest of this iteration.
        HepMC3::GenEvent ge(HepMC3::Units::GEV, HepMC3::Units::MM);
        bool event_retrieved = true;
        #pragma omp critical (reading_HepMCEvent)
        {
          event_retrieved = HepMC3io.read_event(ge);

          // FIXME This is a temp solution to ensure that the event reading
          //       stops when there are no more events in the HepMC file.
          //       Remove this once bugfix is implemented in HepMC.
          if ((ge.particles().size() == 0) && (ge.vertices().size() == 0)) event_retrieved = false;
        }
        if (not event_retrieved) Loop::halt();

        // Translate to HEPUtils event
        get_HEPUtils_event(ge, result);
 
      }

    }
  }

}

#endif
