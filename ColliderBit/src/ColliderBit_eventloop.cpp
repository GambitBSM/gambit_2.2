//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  Functions of ColliderBit event loop.
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
///
///  \author Anders Kvellestad
///          (anders.kvellestad@fys.uio.no)
///  \date   2017 March
///  \date   2018 Jan
///  \date   2018 May
///
///  *********************************************

#include "gambit/Elements/gambit_module_headers.hpp"
#include "gambit/ColliderBit/ColliderBit_eventloop.hpp"

// #define COLLIDERBIT_DEBUG
#define DEBUG_PREFIX "DEBUG: OMP thread " << omp_get_thread_num() << ":  "

namespace Gambit
{

  namespace ColliderBit
  {


    /// LHC Loop Manager
    void operateLHCLoop(MCLoopInfo& result)
    {
      using namespace Pipes::operateLHCLoop;
      static std::streambuf *coutbuf = std::cout.rdbuf(); // save cout buffer for running the loop quietly

      #ifdef COLLIDERBIT_DEBUG
      cout << DEBUG_PREFIX << endl;
      cout << DEBUG_PREFIX << "~~~~ New point! ~~~~" << endl;
      #endif

      // Retrieve run options from the YAML file (or standalone code)
      static bool first = true;
      static bool silenceLoop;
      static std::map<str,int> min_nEvents;
      static std::map<str,int> max_nEvents;
      static std::map<str,int> stoppingres;
      if (first)
      {
        // Should we silence stdout during the loop?
        silenceLoop = runOptions->getValueOrDef<bool>(true, "silenceLoop");

        // Retrieve all the names of all entries in the yaml options node.
        std::vector<str> vec = runOptions->getNames();
        // Step though the names, and accept only those with a "min_nEvents" sub-entry as colliders.
        for (str& name : vec)
        {
          YAML::Node node = runOptions->getNode(name);
          if (not node.IsScalar() and node["min_nEvents"]) result.collider_names.push_back(name);
        }

        // Retrieve the options for each collider.
        for (auto& collider : result.collider_names)
        {
          Options colOptions(runOptions->getValue<YAML::Node>(collider));
          min_nEvents[collider]                                           = colOptions.getValue<int>("min_nEvents");
          max_nEvents[collider]                                           = colOptions.getValue<int>("max_nEvents");
          result.convergence_options[collider].target_stat                = colOptions.getValue<double>("target_fractional_uncert");
          result.convergence_options[collider].stop_at_sys                = colOptions.getValueOrDef<bool>(true, "halt_when_systematic_dominated");
          result.convergence_options[collider].all_analyses_must_converge = colOptions.getValueOrDef<bool>(false, "all_analyses_must_converge");
          result.convergence_options[collider].all_SR_must_converge       = colOptions.getValueOrDef<bool>(false, "all_SR_must_converge");
          result.maxFailedEvents[collider]                                = colOptions.getValueOrDef<int>(1, "maxFailedEvents");
          result.invalidate_failed_points[collider]                       = colOptions.getValueOrDef<bool>(false, "invalidate_failed_points");
          stoppingres[collider]                                           = colOptions.getValueOrDef<int>(200, "events_between_convergence_checks");
          result.analyses[collider]                                       = colOptions.getValueOrDef<std::vector<str>>(std::vector<str>(), "analyses");
          result.event_count[collider]                                    = 0;
          // Check that the nEvents options given make sense.
          if (min_nEvents.at(collider) > max_nEvents.at(collider))
           ColliderBit_error().raise(LOCAL_INFO,"Option min_nEvents is greater than corresponding max_nEvents for collider "
                                                +collider+". Please correct your YAML file.");
          // Check that the analyses all correspond to actual ColliderBit analyses, and sort them into separate maps for each detector.
          for (str& analysis : result.analyses.at(collider))
          {
            result.detector_analyses[collider][getDetector(analysis)].push_back(analysis);
          }
        }
        first = false;
      }

      // Do the base-level initialisation
      Loop::executeIteration(BASE_INIT);

      // Mute stdout during the loop if requested
      if (silenceLoop) std::cout.rdbuf(0);

      // For every collider requested in the yaml file:
      for (auto& collider : result.collider_names)
      {

        // Reset the event_generation_began and exceeded_maxFailedEvents flags
        result.reset_flags();

        // Update the collider
        result.set_current_collider(collider);

        // Initialise the count of the number of generated events.
        result.current_event_count() = 0;

        #ifdef COLLIDERBIT_DEBUG
        cout << DEBUG_PREFIX << "operateLHCLoop: Current collider is " << collider << "." << endl;
        #endif

        piped_invalid_point.check();
        Loop::reset();
        #ifdef COLLIDERBIT_DEBUG
        cout << DEBUG_PREFIX << "operateLHCLoop: Will execute COLLIDER_INIT" << endl;
        #endif
        Loop::executeIteration(COLLIDER_INIT);

        // Any problem during COLLIDER_INIT step?
        piped_warnings.check(ColliderBit_warning());
        piped_errors.check(ColliderBit_error());

        //
        // OMP parallelized sections begin here
        //
        #ifdef COLLIDERBIT_DEBUG
        cout << DEBUG_PREFIX << "operateLHCLoop: Will execute START_SUBPROCESS" << endl;
        #endif
        int currentEvent = 0;
        #pragma omp parallel
        {
          Loop::executeIteration(START_SUBPROCESS);
        }
        // Any problems during the START_SUBPROCESS step?
        piped_warnings.check(ColliderBit_warning());
        piped_errors.check(ColliderBit_error());

        // Convergence loop
        while(currentEvent < max_nEvents.at(collider) and not *Loop::done)
        {
          int eventCountBetweenConvergenceChecks = 0;

          #ifdef COLLIDERBIT_DEBUG
          cout << DEBUG_PREFIX << "Starting main event loop.  Will do " << stoppingres.at(collider) << " events before testing convergence." << endl;
          #endif

          // Main event loop
          #pragma omp parallel
          {
            while(eventCountBetweenConvergenceChecks < stoppingres.at(collider) and
                  currentEvent < max_nEvents.at(collider) and
                  not *Loop::done and
                  not result.exceeded_maxFailedEvents and
                  not piped_warnings.inquire("exceeded maxFailedEvents") and
                  not piped_errors.inquire()
                  )
            {
              result.event_generation_began = true;
              try
              {
                Loop::executeIteration(currentEvent);
                #pragma omp critical
                {
                  currentEvent++;
                  eventCountBetweenConvergenceChecks++;
                }
              }
              catch (std::domain_error& e)
              {
                cout << "\n   Continuing to the next event...\n\n";
              }
            }
          }
          // Update the flag indicating if there have been warnings raised about exceeding the maximum allowed number of failed events
          result.exceeded_maxFailedEvents = result.exceeded_maxFailedEvents or piped_warnings.inquire("exceeded maxFailedEvents");

          // Any problems during the main event loop?
          piped_warnings.check(ColliderBit_warning());
          piped_errors.check(ColliderBit_error());

          #ifdef COLLIDERBIT_DEBUG
          cout << DEBUG_PREFIX << "Did " << eventCountBetweenConvergenceChecks << " events of " << currentEvent << " simulated so far." << endl;
          #endif

          // Break convergence loop if too many events fail
          if(result.exceeded_maxFailedEvents) break;

          // Don't bother with convergence stuff if we haven't passed the minimum number of events yet
          if (currentEvent >= min_nEvents.at(collider))
          {
            #pragma omp parallel
            {
              Loop::executeIteration(COLLECT_CONVERGENCE_DATA);
            }
            // Any problems during the COLLECT_CONVERGENCE_DATA step?
            piped_warnings.check(ColliderBit_warning());
            piped_errors.check(ColliderBit_error());

            Loop::executeIteration(CHECK_CONVERGENCE);
            // Any problems during the CHECK_CONVERGENCE step?
            piped_warnings.check(ColliderBit_warning());
            piped_errors.check(ColliderBit_error());
          }

        }

        // Store the number of generated events
        result.current_event_count() = currentEvent;

        // Skip to next collider if too many events fail
        if(result.exceeded_maxFailedEvents) continue;

        #pragma omp parallel
        {
          Loop::executeIteration(END_SUBPROCESS);
        }
        // Any problems during the END_SUBPROCESS step?
        piped_warnings.check(ColliderBit_warning());
        piped_errors.check(ColliderBit_error());

        //
        // OMP parallelized sections end here
        //

        Loop::executeIteration(COLLIDER_FINALIZE);
      }

      // Nicely thank the loop for being quiet, and restore everyone's vocal chords
      if (silenceLoop) std::cout.rdbuf(coutbuf);

      // Check for exceptions
      piped_invalid_point.check();

      Loop::executeIteration(BASE_FINALIZE);
    }


    /// Store some information about the event generation
    void getLHCEventLoopInfo(map_str_dbl& result)
    {
      using namespace Pipes::getLHCEventLoopInfo;
      result.clear();
      result["did_event_generation"] = double(Dep::RunMC->event_generation_began);
      result["too_many_failed_events"] = double(Dep::RunMC->exceeded_maxFailedEvents);
      for (auto& name : Dep::RunMC->collider_names)
      {
        result["event_count_" + name] = Dep::RunMC->event_count.at(name);
      }
    }


    /// Loop over all analyses and collect them in one place
    void CollectAnalyses(AnalysisDataPointers& result)
    {
      using namespace Pipes::CollectAnalyses;
      static bool first = true;

      // Start with an empty vector
      result.clear();

      #ifdef COLLIDERBIT_DEBUG
        cout << DEBUG_PREFIX << "CollectAnalyses: Dep::ATLASAnalysisNumbers->size()    = " << Dep::ATLASAnalysisNumbers->size() << endl;
        cout << DEBUG_PREFIX << "CollectAnalyses: Dep::CMSAnalysisNumbers->size()      = " << Dep::CMSAnalysisNumbers->size() << endl;
        cout << DEBUG_PREFIX << "CollectAnalyses: Dep::IdentityAnalysisNumbers->size() = " << Dep::IdentityAnalysisNumbers->size() << endl;
      #endif

      // Add results
      if (Dep::ATLASAnalysisNumbers->size() != 0) result.insert(result.end(), Dep::ATLASAnalysisNumbers->begin(), Dep::ATLASAnalysisNumbers->end());
      if (Dep::CMSAnalysisNumbers->size() != 0) result.insert(result.end(), Dep::CMSAnalysisNumbers->begin(), Dep::CMSAnalysisNumbers->end());
      if (Dep::IdentityAnalysisNumbers->size() != 0) result.insert(result.end(), Dep::IdentityAnalysisNumbers->begin(), Dep::IdentityAnalysisNumbers->end());

      // When first called, check that all analyses contain at least one signal region.
      if (first)
      {
        // Loop over all AnalysisData pointers
        for (auto& adp : result)
        {
          if (adp->size() == 0)
          {
            str errmsg;
            errmsg = "The analysis " + adp->analysis_name + " has no signal regions.";
            ColliderBit_error().raise(LOCAL_INFO, errmsg);
          }
        }
        first = false;
      }


      // #ifdef COLLIDERBIT_DEBUG
      // cout << DEBUG_PREFIX << "CollectAnalyses: Current size of 'result': " << result.size() << endl;
      // if (result.size() > 0)
      // {
      //   cout << DEBUG_PREFIX << "CollectAnalyses: Will loop through 'result'..." << endl;
      //   for (auto& adp : result)
      //   {
      //     cout << DEBUG_PREFIX << "CollectAnalyses: 'result' contains AnalysisData pointer to " << adp << endl;
      //     cout << DEBUG_PREFIX << "CollectAnalyses: -- Will now loop over all signal regions in " << adp << endl;
      //     for (auto& sr : adp->srdata)
      //     {
      //       cout << DEBUG_PREFIX << "CollectAnalyses: -- " << adp << " contains signal region: " << sr.sr_label << ", n_signal = " << sr.n_signal << ", n_signal_at_lumi = " << n_signal_at_lumi << endl;
      //     }
      //     cout << DEBUG_PREFIX << "CollectAnalyses: -- Done looping over signal regions in " << adp << endl;
      //   }
      //   cout << DEBUG_PREFIX << "CollectAnalyses: ...Done looping through 'result'." << endl;
      // }
      // #endif
    }


  }

}
