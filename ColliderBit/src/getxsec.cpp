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

      // Read options from yaml file
      static bool first = true;

      const static Finteger inlo = runOptions->getValueOrDef<Finteger>(1, "inlo");                 // specify LO only[0] or complete NLO (slower)[1]
      const static Finteger isq_ng_in = runOptions->getValueOrDef<Finteger>(1, "isq_ng_in");       // specify degenerate [0] or free [1] squark masses
      const static Finteger icoll_in = runOptions->getValueOrDef<Finteger>(1, "icoll_in");         // collider : tevatron[0], lhc[1]
      const static Fdouble energy_in = runOptions->getValueOrDef<Fdouble>(13000.0, "energy_in");  // collider energy in GeV
      const static Finteger i_error_in = runOptions->getValueOrDef<Finteger>(0, "i_error_in");     // with central scale [0] or scale variation [1]

      // ng
      const static bool ng_all = runOptions->getValueOrDef<bool>(false, "ng_all");
      static std::vector<int> ng_ipart1 = runOptions->getValueOrDef<std::vector<int>>(std::vector<int>(), "ng_ipart1");

      // nn
      const static bool nn_all = runOptions->getValueOrDef<bool>(false, "nn_all");
      static std::vector<std::pair<int,int> > nn_ipart1_ipart2 = runOptions->getValueOrDef<std::vector<std::pair<int,int>>>(std::vector<std::pair<int,int>>(), "nn_ipart1_ipart2");

      // ns
      const static bool ns_all = runOptions->getValueOrDef<bool>(false, "ns_all");
      static std::vector<std::pair<int,int> > ns_ipart1_isquark1 = runOptions->getValueOrDef<std::vector<std::pair<int,int>>>(std::vector<std::pair<int,int>>(), "ns_ipart1_isquark1");

      // ll
      const static bool ll_all = runOptions->getValueOrDef<bool>(false, "ll_all");
      static std::vector<int> ll_ipart1 = runOptions->getValueOrDef<std::vector<int>>(std::vector<int>(), "ll_ipart1");

      // tb
      const static bool tb_all = runOptions->getValueOrDef<bool>(false, "tb_all");
      static std::vector<int> tb_ipart1 = runOptions->getValueOrDef<std::vector<int>>(std::vector<int>(), "tb_ipart1");

      // bb
      const static bool bb_all = runOptions->getValueOrDef<bool>(false, "bb_all");
      static std::vector<int> bb_ipart1 = runOptions->getValueOrDef<std::vector<int>>(std::vector<int>(), "bb_ipart1");


      // Use complete process sets?
      if (first) {

        first = false;

        // ng
        if (ng_all) {
          std::vector<int> ng_all_processes = {1, 2, 3, 4, 5, 6, 7, 8};
          ng_ipart1.clear();
          ng_ipart1 = ng_all_processes;
        }

        // nn
        if (nn_all) {
          std::vector<std::pair<int, int> > nn_all_processes = {
            {1,1}, {1,2}, {1,3}, {1,4}, {1,5}, {1,6}, {1,7}, {1,8},
            {2,2}, {2,3}, {2,4}, {2,5}, {2,6}, {2,7}, {2,8},
            {3,3}, {3,4}, {3,5}, {3,6}, {3,7}, {3,8},
            {4,4}, {4,5}, {4,6}, {4,7}, {4,8},
            {5,5}, {5,6}, {5,7}, {5,8},
            {6,6}, {6,7}, {6,8},
            {7,7}, {7,8},
            {8,8}
          };
          nn_ipart1_ipart2.clear();
          nn_ipart1_ipart2 = nn_all_processes;
        }

        // ns
        if (ns_all) {
          std::vector<std::pair<int, int> > ns_all_processes = {
            {1,-5}, {1,-4}, {1,-3}, {1,-2}, {1,-1}, {1,1}, {1,2}, {1,3}, {1,4}, {1,5},
            {2,-5}, {2,-4}, {2,-3}, {2,-2}, {2,-1}, {2,1}, {2,2}, {2,3}, {2,4}, {2,5},
            {3,-5}, {3,-4}, {3,-3}, {3,-2}, {3,-1}, {3,1}, {3,2}, {3,3}, {3,4}, {3,5},
            {4,-5}, {4,-4}, {4,-3}, {4,-2}, {4,-1}, {4,1}, {4,2}, {4,3}, {4,4}, {4,5},
            {5,-5}, {5,-4}, {5,-3}, {5,-2}, {5,-1}, {5,1}, {5,2}, {5,3}, {5,4}, {5,5},
            {6,-5}, {6,-4}, {6,-3}, {6,-2}, {6,-1}, {6,1}, {6,2}, {6,3}, {6,4}, {6,5},
            {7,-5}, {7,-4}, {7,-3}, {7,-2}, {7,-1}, {7,1}, {7,2}, {7,3}, {7,4}, {7,5},
            {8,-5}, {8,-4}, {8,-3}, {8,-2}, {8,-1}, {8,1}, {8,2}, {8,3}, {8,4}, {8,5}
          };
          ns_ipart1_isquark1.clear();
          ns_ipart1_isquark1 = ns_all_processes;
        }

        // ll
        if (ll_all) {
          std::vector<int> ll_all_processes = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13};
          ll_ipart1.clear();
          ll_ipart1 = ll_all_processes;
        }

        // tb
        if (tb_all) {
          std::vector<int> tb_all_processes = {1, 2};
          tb_ipart1.clear();
          tb_ipart1 = tb_all_processes;
        }

        // bb
        if (bb_all) {
          std::vector<int> bb_all_processes = {1, 2};
          bb_ipart1.clear();
          bb_ipart1 = bb_all_processes;
        }


      }

      // Default values
      Finteger ipart1_in = 1;
      Finteger ipart2_in = 2;
      Finteger isquark1_in = 0;
      Finteger isquark2_in = 0;

      // 
      // FIXME: make process loop 
      // 

      // std::string final_state_in = "nn";
      // ipart1_in = nn_ipart1_ipart2.at(0).first;
      // ipart2_in = nn_ipart1_ipart2.at(0).second;

      // std::string final_state_in = "ng";
      // ipart1_in = ng_ipart1.at(0);

      std::string final_state_in = "tb";
      ipart1_in = ng_ipart1.at(0);



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

        ps.inlo = inlo;
        ps.isq_ng_in = isq_ng_in;
        ps.icoll_in = icoll_in;
        ps.energy_in = energy_in;
        ps.i_error_in = i_error_in;
        // ps.final_state_in = final_state_in.c_str();
        ps.final_state_in = final_state_in;
        ps.ipart1_in = ipart1_in;
        ps.ipart2_in = ipart2_in;
        ps.isquark1_in = isquark1_in;
        ps.isquark2_in = isquark2_in;

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
