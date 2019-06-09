//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  Module functions for computing cross-sections
///
///  *********************************************
///
///  Authors (add name and date if you modify):
///
///  \author Pat Scott
///          (p.scott@imperial.ac.uk)
///  \date 2019 Feb, May
///
///  *********************************************

#include "gambit/ColliderBit/ColliderBit_eventloop.hpp"

//#define COLLIDERBIT_DEBUG

namespace Gambit
{

  namespace ColliderBit
  {

    /// Compute a cross-section from Monte Carlo
    void getMCxsec(xsec& result)
    {
      using namespace Pipes::getMCxsec;

      // Don't bother if there are no analyses that will use this.
      if (Dep::RunMC->analyses.empty()) return;

      // Reset the xsec objects on all threads
      if (*Loop::iteration == START_SUBPROCESS)
      {
        result.reset();
      }

      // If we are in the main event loop, count the event towards cross-section normalisation on this thread
      if (*Loop::iteration >= 0)
      {
        result.log_event();
      }

      // Extract the xsecs from the MC on each thread
      if (*Loop::iteration == END_SUBPROCESS && Dep::RunMC->event_generation_began)
      {
        if (not Dep::RunMC->exceeded_maxFailedEvents)
        {
          const double xs_fb = (*Dep::HardScatteringSim)->xsec_pb() * 1000.;
          const double xserr_fb = (*Dep::HardScatteringSim)->xsecErr_pb() * 1000.;
          result.set_xsec(xs_fb, xserr_fb);
          #ifdef COLLIDERBIT_DEBUG
            cout << debug_prefix() << "xs_fb = " << xs_fb << " +/- " << xserr_fb << endl;
          #endif
        }
      }

      // Gather the xsecs from all threads into one
      if (*Loop::iteration == COLLIDER_FINALIZE)
      {
        result.gather_xsecs();
      }

    }

    /// Get a cross-section from NLL-FAST
    void getNLLFastxsec(xsec& result)
    {
      using namespace Pipes::getNLLFastxsec;

      // Don't bother if there are no analyses that will use this.
      if (Dep::RunMC->analyses.empty()) return;

      // Reset the xsec objects on all threads
      if (*Loop::iteration == START_SUBPROCESS)
      {
        result.reset();
      }

      // If we are in the main event loop, count the event towards cross-section normalisation on this thread
      if (*Loop::iteration >= 0)
      {
        result.log_event();
      }

      // Set the xsec and its error, and gather event counts from all threads.
      if (*Loop::iteration == COLLIDER_FINALIZE)
      {
        double xs_fb = 0.1;             // replace with xsec from NLL-Fast
        double xserr_fb = 0.1 * xs_fb;  // or whatever
        result.set_xsec(xs_fb, xserr_fb);
        result.gather_num_events();
      }
    }

    /// A function that reads the total cross-section from the input file, but builds up the number of events from the event loop
    void getYAMLxsec(xsec& result)
    {
      using namespace Pipes::getYAMLxsec;

      // Don't bother if there are no analyses that will use this.
      if (Dep::RunMC->analyses.empty()) return;

      // Retrieve the total cross-section and cross-section error
      const static double xsec_fb = 1000 * runOptions->getValue<double>("xsec_pb");
      const static double xsec_fractional_uncert = runOptions->getValue<double>("xsec_fractional_uncert");

      // Reset the xsec objects on all threads
      if (*Loop::iteration == START_SUBPROCESS)
      {
        result.reset();
      }

      // If we are in the main event loop, count the event towards cross-section normalisation on this thread
      if (*Loop::iteration >= 0)
      {
        result.log_event();
      }

      // Set the xsec and its error
      if (*Loop::iteration == COLLIDER_FINALIZE)
      {
        result.set_xsec(xsec_fb, xsec_fractional_uncert*xsec_fb);
        result.gather_num_events();
      }

    }


    /// Get a cross-section from Prospino
    void getProspinoxsec(xsec& result)
    {
      using namespace Pipes::getProspinoxsec;

      // // Don't bother if there are no analyses that will use this.
      // if (Dep::RunMC->analyses.empty()) return;

      cout << "DEBUG: getProspinoxsec: Loop iteration:" << *Loop::iteration << endl;

      // Reset the xsec object on the main thread (other threads do not matter)
      if (*Loop::iteration == BASE_INIT)
      {
        result.reset();

        // Testing...
        cout << "DEBUG: getProspinoxsec: Requesting dependency BE::prospino_LHC_xsec..." << endl;


        // Get an SLHA1 object for Prospino.
        const SLHAstruct& slha = Dep::MSSM_spectrum->getSLHAea(1);

        // Get the GAMBIT model parameters for Prospino
        const param_map_type& model_params = Param;

        // Create struct with settings Prospino
        prospino_settings ps;

        ps.inlo = 0;
        ps.isq_ng_in = 1;
        ps.icoll_in = 1;
        ps.energy_in = 13000.;
        ps.i_error_in = 0;
        ps.final_state_in = "nn";
        ps.ipart1_in = 1;
        ps.ipart2_in = 2;
        ps.isquark1_in = 0;
        ps.isquark2_in = 0;

        // inlo = runOptions->getValueOrDef<Finteger>(1, "inlo");                 // specify LO only[0] or complete NLO (slower)[1]
        // isq_ng_in = runOptions->getValueOrDef<Finteger>(1, "isq_ng_in");       // specify degenerate [0] or free [1] squark masses
        // icoll_in = runOptions->getValueOrDef<Finteger>(1, "icoll_in");         // collider : tevatron[0], lhc[1]
        // energy_in = runOptions->getValueOrDef<Fdouble>(13000.0, "energy_in");  // collider energy in GeV
        // i_error_in = runOptions->getValueOrDef<Finteger>(0, "i_error_in");     // with central scale [0] or scale variation [1]
        
        // final_state_in = runOptions->getValueOrDef<std::string>("nn", "final_state_in"); // select process
        // ipart1_in = runOptions->getValueOrDef<Finteger>(1, "ipart1_in");      //
        // ipart2_in = runOptions->getValueOrDef<Finteger>(2, "ipart2_in");      //
        // isquark1_in = runOptions->getValueOrDef<Finteger>(0, "isquark1_in");  //
        // isquark2_in = runOptions->getValueOrDef<Finteger>(0, "isquark2_in");  //


        // Call Prospino and get the result in a map<string,double>
        map_str_dbl prospino_output = BEreq::prospino_LHC_xsec(slha, model_params, ps);

        cout << "DEBUG: getProspinoxsec: got this result from Prospino:" << endl;

        for (auto& element : prospino_output)
        {
          cout << "DEBUG: getProspinoxsec: " << element.first << " = " << element.second << endl;
        }

        // set result
        result.set_xsec(prospino_output.at("NLO_ms[pb]"), prospino_output.at("NLO_rel_error"));
      }

      // If we are in the main event loop, count the event towards cross-section normalisation on this thread
      if (*Loop::iteration >= 0)
      {
        result.log_event();
      }

      // Set the xsec and its error
      if (*Loop::iteration == COLLIDER_FINALIZE)
      {
        result.gather_num_events();
      }


    }



  }
}
