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

// #define COLLIDERBIT_DEBUG
#define DEBUG_PREFIX "DEBUG: OMP thread " << omp_get_thread_num() << ":  "

namespace Gambit
{

  namespace ColliderBit
  {

    // ======= Utility functions =======


    /// Helper function that takes a cross-section value in fb or pb, 
    /// along with an absolute or relative uncertainty, and returns the 
    /// xsec and absolute uncertainty in fb.
    std::pair<double,double> convert_xsecs_to_fb(double input_xsec, double input_xsec_uncert, str input_unit, bool input_fractional_uncert)
    {
      double xsec_fb;
      double xsec_uncert_fb;

      if (input_unit == "fb" && !input_fractional_uncert)
      {
        xsec_fb = input_xsec;
        xsec_uncert_fb = input_xsec_uncert;
      }
      else if (input_unit == "fb" && input_fractional_uncert)
      {
        xsec_fb = input_xsec;
        xsec_uncert_fb = input_xsec_uncert * xsec_fb;
      }
      else if (input_unit == "pb" && !input_fractional_uncert)
      {
        xsec_fb = input_xsec * 1000.;
        xsec_uncert_fb = input_xsec_uncert * 1000.;
      }
      else if (input_unit == "pb" && input_fractional_uncert)
      {
        xsec_fb = input_xsec * 1000.;
        xsec_uncert_fb = input_xsec_uncert * xsec_fb;
      }
      else
      {
        ColliderBit_error().raise(LOCAL_INFO, "Unknown combination of options for function convert_xsecs_to_fb.");
      }      

      cout << "DEBUG: returning xsec [fb] = " << xsec_fb << " +/- " << xsec_uncert_fb << endl;
      return std::make_pair(xsec_fb, xsec_uncert_fb);
    }


    // ======= Module functions =======


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
            cout << DEBUG_PREFIX << "xs_fb = " << xs_fb << " +/- " << xserr_fb << endl;
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

      static std::pair<str,str> xsec_pnames;
      static str input_unit; 
      static bool input_fractional_uncert = false;


      static bool first = true;
      if (*Loop::iteration == BASE_INIT)
      {

        if (first)
        {
          // Determine the correct combination of parameters
          if ((runOptions->hasKey("xsec_fb")) && (runOptions->hasKey("xsec_uncert_fb")))
          {
            xsec_pnames.first = "xsec_fb";
            xsec_pnames.second = "xsec_uncert_fb";
            input_unit = "fb";
            input_fractional_uncert = false;
          }
          else if ((runOptions->hasKey("xsec_fb")) && (runOptions->hasKey("xsec_fractional_uncert")))
          {
            xsec_pnames.first = "xsec_fb";
            xsec_pnames.second = "xsec_fractional_uncert";
            input_unit = "fb";
            input_fractional_uncert = true;
          }
          else if ((runOptions->hasKey("xsec_pb")) && (runOptions->hasKey("xsec_uncert_pb")))
          {
            xsec_pnames.first = "xsec_pb";
            xsec_pnames.second = "xsec_uncert_pb";
            input_unit = "pb";
            input_fractional_uncert = false;
          }
          else if ((runOptions->hasKey("xsec_pb")) && (runOptions->hasKey("xsec_fractional_uncert")))
          {
            xsec_pnames.first = "xsec_pb";
            xsec_pnames.second = "xsec_fractional_uncert";
            input_unit = "pb";
            input_fractional_uncert = true;
          }
          else
          {
            std::stringstream errmsg_ss;
            errmsg_ss << "Unknown combination of options for function getYAMLxsec." << endl;
            errmsg_ss << "Needs one of the following sets of option names:" << endl;
            errmsg_ss << "  xsec_fb, xsec_uncert_fb" << endl;
            errmsg_ss << "  xsec_fb, xsec_fractional_uncert" << endl;
            errmsg_ss << "  xsec_pb, xsec_uncert_pb" << endl;
            errmsg_ss << "  xsec_pb, xsec_fractional_uncert" << endl;
            ColliderBit_error().raise(LOCAL_INFO, errmsg_ss.str());
          }

          first = false;
        }
      }


      // Retrieve the total cross-section and cross-section error
      const static double input_xsec = runOptions->getValue<double>(xsec_pnames.first);
      const static double input_xsec_uncert = runOptions->getValue<double>(xsec_pnames.second);

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
        std::pair<double,double> temp = convert_xsecs_to_fb(input_xsec, input_xsec_uncert, input_unit, input_fractional_uncert);
        double xsec_fb = temp.first;
        double xsec_uncert_fb = temp.second;

        result.set_xsec(xsec_fb, xsec_uncert_fb);
        result.gather_num_events();
      }

    }


    /// A function that reads a list of (SLHA file, total cross-section) pairs from the input YAML file
    void getYAMLxsec_SLHA(xsec& result)
    {
      using namespace Pipes::getYAMLxsec_SLHA;

      // Don't bother if there are no analyses that will use this.
      if (Dep::RunMC->analyses.empty()) return;

      static std::pair<str,str> xsec_pnames;
      static str input_unit; 
      static bool input_fractional_uncert = false;

      static bool first = true;
      if (*Loop::iteration == BASE_INIT)
      {

        if (first)
        {
          // Determine the correct combination of parameters
          if ((runOptions->hasKey("xsec_fb")) && (runOptions->hasKey("xsec_uncert_fb")))
          {
            xsec_pnames.first = "xsec_fb";
            xsec_pnames.second = "xsec_uncert_fb";
            input_unit = "fb";
            input_fractional_uncert = false;
          }
          else if ((runOptions->hasKey("xsec_fb")) && (runOptions->hasKey("xsec_fractional_uncert")))
          {
            xsec_pnames.first = "xsec_fb";
            xsec_pnames.second = "xsec_fractional_uncert";
            input_unit = "fb";
            input_fractional_uncert = true;
          }
          else if ((runOptions->hasKey("xsec_pb")) && (runOptions->hasKey("xsec_uncert_pb")))
          {
            xsec_pnames.first = "xsec_pb";
            xsec_pnames.second = "xsec_uncert_pb";
            input_unit = "pb";
            input_fractional_uncert = false;
          }
          else if ((runOptions->hasKey("xsec_pb")) && (runOptions->hasKey("xsec_fractional_uncert")))
          {
            xsec_pnames.first = "xsec_pb";
            xsec_pnames.second = "xsec_fractional_uncert";
            input_unit = "pb";
            input_fractional_uncert = true;
          }
          else
          {
            std::stringstream errmsg_ss;
            errmsg_ss << "Unknown combination of options for function getYAMLxsec_SLHA." << endl;
            errmsg_ss << "Needs one of the following sets of option names:" << endl;
            errmsg_ss << "  xsec_fb, xsec_uncert_fb" << endl;
            errmsg_ss << "  xsec_fb, xsec_fractional_uncert" << endl;
            errmsg_ss << "  xsec_pb, xsec_uncert_pb" << endl;
            errmsg_ss << "  xsec_pb, xsec_fractional_uncert" << endl;
            ColliderBit_error().raise(LOCAL_INFO, errmsg_ss.str());
          }

          first = false;
        }
      }

      // Get the filename, look for it in the xsec and uncertainty lists
      const static YAML::Node colNode_xsec = runOptions->getValue<YAML::Node>(xsec_pnames.first);
      const static Options colOptions_xsec(colNode_xsec);
      const static YAML::Node colNode_uncert = runOptions->getValue<YAML::Node>(xsec_pnames.second);
      const static Options colOptions_uncert(colNode_uncert);
      static str filename;

      if (*Loop::iteration == BASE_INIT)
      {
        // Update the SLHA filename
        filename = Dep::SLHAFileNameAndContent->first;

        // Look for the filename in the xsec lists
        if (!colOptions_xsec.hasKey(filename)) piped_invalid_point.request(str("No cross-section found for SLHA file ").append(filename));
        if (!colOptions_uncert.hasKey(filename)) piped_invalid_point.request(str("No fractional cross-section uncertainty found for SLHA file ").append(filename));
      }

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
        double input_xsec = colOptions_xsec.getValue<double>(filename);
        double input_xsec_uncert = colOptions_uncert.getValue<double>(filename);

        std::pair<double,double> temp = convert_xsecs_to_fb(input_xsec, input_xsec_uncert, input_unit, input_fractional_uncert);
        double xsec_fb = temp.first;
        double xsec_uncert_fb = temp.second;

        result.set_xsec(xsec_fb, xsec_uncert_fb);
        result.gather_num_events();
      }

    }  // end getYAMLxsec_SLHA



    /// A function that assigns a total cross-sections directly from the scan parameters
    /// (for models CB_SLHA_simpmod_scan_model and CB_SLHA_scan_model)
    void getYAMLxsec_param(xsec& result)
    {
      using namespace Pipes::getYAMLxsec_param;

      // Don't bother if there are no analyses that will use this.
      if (Dep::RunMC->analyses.empty()) return;

      static std::vector<str> pnames;
      static std::pair<str,str> xsec_pnames;

      static str input_unit; 
      static bool input_fractional_uncert = false;

      static bool first = true;
      if (*Loop::iteration == BASE_INIT)
      {

        if (first)
        {

          // Get all parameter names
          for (const auto& parname_parptr_pair : Param)
          {
            pnames.push_back(parname_parptr_pair.first);
          }

          // Determine the correct combination of parameters
          if ((std::find(pnames.begin(), pnames.end(), "xsec_fb") != pnames.end()) 
               && (std::find(pnames.begin(), pnames.end(), "xsec_uncert_fb") != pnames.end()))
          {
            xsec_pnames.first = "xsec_fb";
            xsec_pnames.second = "xsec_uncert_fb";
            input_unit = "fb";
            input_fractional_uncert = false;
          }
          else if ((std::find(pnames.begin(), pnames.end(), "xsec_fb") != pnames.end()) 
                    && (std::find(pnames.begin(), pnames.end(), "xsec_fractional_uncert") != pnames.end()))
          {
            xsec_pnames.first = "xsec_fb";
            xsec_pnames.second = "xsec_fractional_uncert";
            input_unit = "fb";
            input_fractional_uncert = true;
          }
          else if ((std::find(pnames.begin(), pnames.end(), "xsec_pb") != pnames.end()) 
                    && (std::find(pnames.begin(), pnames.end(), "xsec_uncert_pb") != pnames.end()))
          {
            xsec_pnames.first = "xsec_pb";
            xsec_pnames.second = "xsec_uncert_pb";
            input_unit = "pb";
            input_fractional_uncert = false;
          }
          else if ((std::find(pnames.begin(), pnames.end(), "xsec_pb") != pnames.end()) 
                    && (std::find(pnames.begin(), pnames.end(), "xsec_fractional_uncert") != pnames.end()))
          {
            xsec_pnames.first = "xsec_pb";
            xsec_pnames.second = "xsec_fractional_uncert";
            input_unit = "pb";
            input_fractional_uncert = true;
          }
          else
          {
            std::stringstream errmsg_ss;
            errmsg_ss << "Unknown combination of parameters for function getYAMLxsec_param." << endl;
            errmsg_ss << "Needs one of the following sets of parameter names:" << endl;
            errmsg_ss << "  xsec_fb, xsec_uncert_fb" << endl;
            errmsg_ss << "  xsec_fb, xsec_fractional_uncert" << endl;
            errmsg_ss << "  xsec_pb, xsec_uncert_pb" << endl;
            errmsg_ss << "  xsec_pb, xsec_fractional_uncert" << endl;
            ColliderBit_error().raise(LOCAL_INFO, errmsg_ss.str());
          }

          first = false;
        }
      }

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
        double input_xsec = *Param.at(xsec_pnames.first);
        double input_xsec_uncert = *Param.at(xsec_pnames.second); 

        std::pair<double,double> temp = convert_xsecs_to_fb(input_xsec, input_xsec_uncert, input_unit, input_fractional_uncert);
        double xsec_fb = temp.first;
        double xsec_uncert_fb = temp.second;

        result.set_xsec(xsec_fb, xsec_uncert_fb);
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
      Finteger ipart2_in = 1;
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

      // std::string final_state_in = "tb";
      // ipart1_in = ng_ipart1.at(0);

      // Reset the xsec object on the main thread (other threads do not matter)
      if (*Loop::iteration == BASE_INIT)
      {
        result.reset();

        // Testing...
        cout << "DEBUG: getProspinoxsec: Requesting dependency BE::prospino_LHC_xsec..." << endl;

        // Accumulators
        double combined_xsec_pb = 0.0;
        double combined_xsec_err = 0.0;
        double combined_xsec_rel_err = 0.0;

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

        // ng
        cout << "DEBUG: getProspinoxsec: " << endl;
        cout << "DEBUG: getProspinoxsec: --------------------------------" << endl;
        cout << "DEBUG: getProspinoxsec: Process: ng" << endl;
        cout << "DEBUG: getProspinoxsec: --------------------------------" << endl;
        cout << "DEBUG: getProspinoxsec: " << endl;
        for (int& ip1 : ng_ipart1) {

          ps.final_state_in = "ng";
          ps.ipart1_in = ip1;
          ps.ipart2_in = 1;
          ps.isquark1_in = 0;
          ps.isquark2_in = 0;

          // Call Prospino and get the result in a map<string,double>
          map_str_dbl prospino_output = BEreq::prospino_LHC_xsec(slha, model_params, ps);

          cout << "DEBUG: getProspinoxsec: process: ng" << endl;
          cout << "DEBUG: getProspinoxsec: ipart1, ipart2: " << ps.ipart1_in << ", " << ps.ipart2_in << endl;
          cout << "DEBUG: getProspinoxsec: isquark1, isquark2: " << ps.isquark1_in << ", " << ps.isquark2_in << endl;
          cout << "DEBUG: getProspinoxsec: got this result from Prospino:" << endl;
          for (auto& element : prospino_output)
          {
            cout << "DEBUG: getProspinoxsec: " << element.first << " = " << element.second << endl;
          }

          double xs_pb = prospino_output.at("NLO_ms[pb]");
          double xs_rel_err = prospino_output.at("NLO_rel_error");
          double xs_err = xs_pb * xs_rel_err;

          combined_xsec_pb += xs_pb;
          combined_xsec_err = sqrt(pow(combined_xsec_err,2) + pow(xs_err,2));
          combined_xsec_rel_err = combined_xsec_err / combined_xsec_pb;

          result.set_xsec(combined_xsec_pb, combined_xsec_rel_err);
        }


        // nn
        cout << "DEBUG: getProspinoxsec: " << endl;
        cout << "DEBUG: getProspinoxsec: --------------------------------" << endl;
        cout << "DEBUG: getProspinoxsec: Process: nn" << endl;
        cout << "DEBUG: getProspinoxsec: --------------------------------" << endl;
        cout << "DEBUG: getProspinoxsec: " << endl;
        for (std::pair<int,int>& ip1_ip2_pair : nn_ipart1_ipart2) {

          ps.final_state_in = "nn";
          ps.ipart1_in = ip1_ip2_pair.first;
          ps.ipart2_in = ip1_ip2_pair.second;
          ps.isquark1_in = 0;
          ps.isquark2_in = 0;

          // Call Prospino and get the result in a map<string,double>
          map_str_dbl prospino_output = BEreq::prospino_LHC_xsec(slha, model_params, ps);

          cout << "DEBUG: getProspinoxsec: process: nn" << endl;
          cout << "DEBUG: getProspinoxsec: ipart1, ipart2: " << ps.ipart1_in << ", " << ps.ipart2_in << endl;
          cout << "DEBUG: getProspinoxsec: isquark1, isquark2: " << ps.isquark1_in << ", " << ps.isquark2_in << endl;
          cout << "DEBUG: getProspinoxsec: got this result from Prospino:" << endl;
          for (auto& element : prospino_output)
          {
            cout << "DEBUG: getProspinoxsec: " << element.first << " = " << element.second << endl;
          }

          double xs_pb = prospino_output.at("NLO_ms[pb]");
          double xs_rel_err = prospino_output.at("NLO_rel_error");
          double xs_err = xs_pb * xs_rel_err;

          combined_xsec_pb += xs_pb;
          combined_xsec_err = sqrt(pow(combined_xsec_err,2) + pow(xs_err,2));
          combined_xsec_rel_err = combined_xsec_err / combined_xsec_pb;

          result.set_xsec(combined_xsec_pb, combined_xsec_rel_err);
        }

        // 
        // Add more process blocks here...
        //

      } // end iteration BASE_INIT

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

    } // end getProspinoxsec



    /// Get cross-section info as map_str_dbl (for simple printing)
    void getXsecInfoMap(map_str_dbl& result)
    {
      using namespace Pipes::getXsecInfoMap;

      // @todo Do we need this to ensure that the result map is always of the same length (for the printer)?
      // // Append the xsec info for the current collider to the result map
      // if (*Loop::iteration == COLLIDER_INIT)
      // {
      //   xsec empty_xs;
      //   for(auto s_d_pair : empty_xs.get_content_as_map())
      //   {
      //     std::string new_key(Dep::RunMC->current_collider());
      //     new_key.append("__").append(s_d_pair.first);
      //     result[new_key] = s_d_pair.second;
      //   }
      // }

      // Append the xsec info for the current collider to the result map
      if (*Loop::iteration == COLLIDER_FINALIZE)
      {
        const xsec& xs = (*Dep::CrossSection);
        for(auto s_d_pair : xs.get_content_as_map())
        {
          std::string new_key(Dep::RunMC->current_collider());
          new_key.append("__").append(s_d_pair.first);
          result[new_key] = s_d_pair.second;
        }
      }
    }  // end getXsecInfoMap

  }
}
