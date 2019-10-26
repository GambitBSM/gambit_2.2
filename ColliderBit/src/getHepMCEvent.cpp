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
///  \author Tomek Procter
///           (tsp116@ic.ac.uk)
///  \date 2019 Octobet
///  *********************************************

#include "gambit/cmake/cmake_variables.hpp"

#ifndef EXCLUDE_HEPMC

#include "gambit/ColliderBit/ColliderBit_eventloop.hpp"
#include "gambit/Utils/util_functions.hpp"
#include "HepMC3/ReaderAsciiHepMC2.h"
#include "gambit/ColliderBit/colliders/Pythia8/Py8EventConversions.hpp"

#include "HepMC3/GenEvent.h"
#include "HepMC3/GenParticle.h"


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
      const static str hepmc_filename = runOptions->getValue<str>("hepmc_filename");
      static bool first = true;
      if (first)
      {
        if (not Utils::file_exists(hepmc_filename)) throw std::runtime_error("HepMC event file "+hepmc_filename+" not found.  Quitting...");
        first = false;
      }
      static HepMC3::ReaderAsciiHepMC2 hepmcio(hepmc_filename);

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
        event_retrieved = hepmcio.read_event(ge);

        // FIXME This is a temp solution to ensure that the event reading
        //       stops when there are no more events in the HepMC file.
        //       Remove this once bugfix is implemented in HepMC.
        if ((ge.particles().size() == 0) && (ge.vertices().size() == 0)) event_retrieved = false;
      }
      if (not event_retrieved) Loop::halt();

      // Translate to HEPUtils event

      //A couple of things need to be done before we pass to the event conversion function: set the
      //weight and 
      result.set_weight(ge.weight());

      //annoyingly, I have to add this extra step here, as I need an explicitly const vector of const particles,
      //and I can't generate that without first making a const HepMC3::GenEvent.
      const HepMC3::GenEvent ge2 = ge;
      const std::vector<HepMC3::ConstGenParticlePtr> particles = ge2.particles();

      //Call the unified HEPMC/Pythia event converter:
      //n.b. the 0.4 is the antiktR -> its always been hardcoded in CBS, should it be?
      Gambit::ColliderBit::convertParticleEvent(particles, result, 0.4);
    }

  }

}

#endif
