//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  ColliderBit event loop functions returning
///  collider Monte Carlo events.
///
///  *********************************************
///
///  Authors (add name and date if you modify):
///
///  \author Abram Krislock
///          (a.m.b.krislock@fys.uio.no)
///
///  \author Aldo Saavedra
///
///  \author Andy Buckley
///
///  \author Chris Rogan
///          (crogan@cern.ch)
///  \date 2014 Aug
///  \date 2015 May
///
///  \author Pat Scott
///          (p.scott@imperial.ac.uk)
///  \date 2015 Jul
///  \date 2018 Jan
///  \date 2019 Jan
///  \date 2019 May
///
///  \author Anders Kvellestad
///          (anders.kvellestad@fys.uio.no)
///  \date   2017 March
///  \date   2018 Jan
///  \date   2018 May
///
///  \author Tomas Gonzalo
///          (tomas.gonzalo@monash.edu)
///  \date   2019 June
///
///  *********************************************

#include "gambit/ColliderBit/ColliderBit_eventloop.hpp"
#include "gambit/ColliderBit/colliders/Pythia8/Py8EventConversions.hpp"

#define COLLIDERBIT_DEBUG

namespace Gambit
{

  namespace ColliderBit
  {

    /// Generate a hard scattering event with Pythia
    template<typename PythiaT, typename EventT>
    void generateEventPy8Collider(HEPUtils::Event& event,
                                  const MCLoopInfo& RunMC,
                                  const Py8Collider<PythiaT,EventT>& HardScatteringSim,
                                  const int iteration,
                                  void(*wrapup)(),
                                  const Options& runOptions)
    {
      static int nFailedEvents;
      thread_local EventT pythia_event;

      // If the event loop has not been entered yet, reset the counter for the number of failed events
      if (iteration == BASE_INIT)
      {
        // Although BASE_INIT should never be called in parallel, we do this in omp_critical just in case.
        #pragma omp critical (pythia_event_failure)
        {
          nFailedEvents = 0;
        }
        return;
      }

      // If in any other special iteration, do nothing
      if (iteration < BASE_INIT) return;

      // Reset the Pythia and HEPUtils events
      pythia_event.clear();
      event.clear();

      // Attempt (possibly repeatedly) to generate an event
      while(nFailedEvents <= RunMC.current_maxFailedEvents())
      {
        try
        {
          HardScatteringSim.nextEvent(pythia_event);
          break;
        }
        catch (typename Py8Collider<PythiaT,EventT>::EventGenerationError& e)
        {
          #ifdef COLLIDERBIT_DEBUG
          cout << debug_prefix() << "Py8Collider::EventGenerationError caught in generateEventPy8Collider. Check the ColliderBit log for event details." << endl;
          #endif
          #pragma omp critical (pythia_event_failure)
          {
            // Update global counter
            nFailedEvents += 1;
            // Store Pythia event record in the logs
            std::stringstream ss;
            pythia_event.list(ss, 1);
            logger() << LogTags::debug << "Py8Collider::EventGenerationError error caught in generateEventPy8Collider. Pythia record for event that failed:\n" << ss.str() << EOM;
          }
        }
      }

      // Wrap up event loop if too many events fail.
      if(nFailedEvents > RunMC.current_maxFailedEvents())
      {
        piped_warnings.request(LOCAL_INFO,"exceeded maxFailedEvents");
        wrapup();
        return;
      }

      #ifndef EXCLUDE_HEPMC
        // Print the event to HepMC
        long long maxIterationsForHepMC = runOptions.getValueOrDef<long long>(50000, "max_HepMC_events");;
        if (iteration <= maxIterationsForHepMC)
        {
          if (runOptions.getValueOrDef<bool>(false, "drop_HepMC_file") )
          {
            cout << "printing HepMC" << endl;

            // Interface for conversion from Pythia8::Event to HepMC event.
//            HepMC3::Pythia8ToHepMC ToHepMC;

            // Specify file where HepMC events will be stored.
//            HepMC3::IO_GenEvent ascii_io("GAMBIT_HepMC.dat", std::ios::out);

            // Construct new empty HepMC event and fill it.
            // Units will be as chosen for HepMC build; but can be changed
            // by arguments, e.g. GenEvt( HepMC::Units::GEV, HepMC::Units::MM)
 //           HepMC3::GenEvent* hepmcevt = new HepMC3::GenEvent();
 //           ToHepMC.fill_next_event( pythia, hepmcevt );

            // Write the HepMC event to file. Done with it.
  //          #pragma once
 //           ascii_io << hepmcevt;
//            delete hepmcevt;
          }

        }
      #endif

      // Attempt to convert the Pythia event to a HEPUtils event
      try
      {
        if (HardScatteringSim.partonOnly)
          convertPartonEvent(pythia_event, event, HardScatteringSim.antiktR);
        else
          convertParticleEvent(pythia_event, event, HardScatteringSim.antiktR);
      }
      // No good.
      catch (Gambit::exception& e)
      {
        #ifdef COLLIDERBIT_DEBUG
          cout << debug_prefix() << "Gambit::exception caught during event conversion in generateEventPy8Collider. Check the ColliderBit log for details." << endl;
        #endif

        #pragma omp critical (event_conversion_error)
        {
          // Store Pythia event record in the logs
          std::stringstream ss;
          pythia_event.list(ss, 1);
          logger() << LogTags::debug << "Gambit::exception caught in generateEventPy8Collider. Pythia record for event that failed:\n" << ss.str() << EOM;
        }

        str errmsg = "Bad point: generateEventPy8Collider caught the following runtime error: ";
        errmsg    += e.what();
        piped_invalid_point.request(errmsg);
        wrapup();
        return;
      }

    }

    /// Generate a hard scattering event with a specific Pythia
    #define GET_PYTHIA_EVENT(NAME)                               \
    void NAME(HEPUtils::Event& result)                           \
    {                                                            \
      using namespace Pipes::NAME;                               \
      generateEventPy8Collider(result, *Dep::RunMC,              \
       *Dep::HardScatteringSim, *Loop::iteration, Loop::wrapup,  \
       *runOptions); \
    }

  }

}
