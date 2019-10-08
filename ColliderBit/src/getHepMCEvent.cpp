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
      const static str HepMC_filename = runOptions->getValueOrDef<str>("", "hepmc_filename");
      static int HepMC_file_version = -1;

      static bool first = true;
      if (first)
      {
        if (not Utils::file_exists(HepMC_filename)) throw std::runtime_error("HepMC event file "+HepMC_filename+" not found. Quitting...");

        // Figure out if the file is HepMC2 or HepMC3
        std::ifstream infile(HepMC_filename);
        if (infile.good())
        {
          std::string line;
          while(std::getline(infile, line))
          {
            // Skip blank lines
            if(line == "") continue;

            // We look for "HepMC::Version 2" or "HepMC::Version 3", 
            // so we only need the first 16 characters of the line
            std::string short_line = line.substr(0,16);

            if (short_line == "HepMC::Version 2")
            {
              HepMC_file_version = 2;
              break;
            }
            else if (short_line == "HepMC::Version 3")
            {
              HepMC_file_version = 3;
              break;
            }
            else
            {
              throw std::runtime_error("Could not determine HepMC version from the string '"+short_line+"' extracted from the line '"+line+"'. Quitting...");
            }
          }
        }
        first = false;
      }

      if(HepMC_file_version != 2 and HepMC_file_version != 3)
      {
        throw std::runtime_error("Failed to determine HepMC version for input file "+HepMC_filename+". Quitting...");
      }

      static HepMC3::Reader *HepMCio;

      // Initialize the reader on the first iteration
      if (*Loop::iteration == BASE_INIT)
      {
        if (HepMC_file_version == 2)
          HepMCio = new HepMC3::ReaderAsciiHepMC2(HepMC_filename);
        else
          HepMCio = new HepMC3::ReaderAscii(HepMC_filename);
      }

      // Delete the reader in the last iteration
      if (*Loop::iteration == BASE_FINALIZE)
        delete HepMCio;

      // Don't do anything else during special iterations
      if (*Loop::iteration < 0) return;

      #ifdef COLLIDERBIT_DEBUG
        cout << "Event number: " << *Loop::iteration << endl;
      #endif

      // Attempt to read the next HepMC event as a HEPUtils event. If there are no more events, wrap up the loop and skip the rest of this iteration.
      HepMC3::GenEvent ge(HepMC3::Units::GEV, HepMC3::Units::MM);
      bool event_retrieved = true;
      #pragma omp critical (reading_HepMCEvent)
      {
        event_retrieved = HepMCio->read_event(ge);
     
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

#endif
