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


    /// Dummy function for testing getProcessCrossSections
    xsec dummyXsecFunction(const SLHAstruct& slha, const PID_pair& pids)
    {
      xsec xs_result;

      // Calculate cross section
      double xs_pb = 3.1415;
      double xs_rel_err = 0.01;
      // double xs_err = xs_pb * xs_rel_err;

      // Save result in xs_result
      xs_result.set_xsec(xs_pb, xs_rel_err);

      // Construct info string of the form "PID1:<PID1>, PID2:<PID2>"
      std::stringstream info_ss;
      info_ss << "PID1:" << pids.first << ", " << "PID2:" << pids.second;
      xs_result.set_info_string(info_ss.str());

      return xs_result;
    }


    // ======= Module functions =======


    /// Get a map between Pythia process codes and cross-sections
    void getProcessCrossSections(map_int_xsec& result)
    {
      using namespace Pipes::getProcessCrossSections;

      if(*Loop::iteration == COLLIDER_INIT_OMP)
      {
        result.clear();
      }

      if(*Loop::iteration == XSEC_CALCULATION)
      {
        cout << DEBUG_PREFIX << "getProcessCrossSections: it = XSEC_CALCULATION, ProcessCodes.size() = " << Dep::ProcessCodes->size() << endl;          
        cout << DEBUG_PREFIX << "getProcessCrossSections: it = XSEC_CALCULATION, ProcessCodeToPIDPairsMap.size() = " << Dep::ProcessCodeToPIDPairsMap->size() << endl;          

        // Get an SLHA1 object
        const SLHAstruct& slha = Dep::MSSM_spectrum->getSLHAea(1);

        // Loop over all active processes and construct the cross-section map (result)
        for (size_t i = 0; i != Dep::ProcessCodes->size(); ++i)
        {
          int pcode = Dep::ProcessCodes->at(i);

          // Accumulate all the cross-sections corresponding to the process code pcode
          xsec xs_combined;

          // Get iterator bounds (as a pair) over the multimap entries that match the key pcode
          auto mm_range = Dep::ProcessCodeToPIDPairsMap->equal_range(pcode);

          // Loop over these elements in the multimap
          for (auto mm_it = mm_range.first; mm_it != mm_range.second; ++mm_it)
          {
            PID_pair pids = mm_it->second;

            // Call cross-section calculator
            xsec xs = dummyXsecFunction(slha, pids);
            // @todo Do we need to add a call with reversed PIDs, e.g. (PID2,PID1)? Depends on the calculator!

            // Accumulate result
            xs_combined.sum_xsecs(xs);
          }

          // Construct info string of the form "ProcessCode:<pcode>"
          std::stringstream info_ss;
          info_ss << "ProcessCode:" << pcode;
          xs_combined.set_info_string(info_ss.str());

          // Save combined cross-section for this process code in the result map
          result[pcode] = xs_combined;
        }


      // // Fill NLO cross sections for turned-on processes
      // for(vector<int>::iterator it = procs.begin(); it != procs.end(); ++it){
      //   // Process code
      //   int process = *it;
      //   // Loop over and sum NLO cross sections from all PID combinations belonging to this process number
      //   double NLOxsec = 0;
      //   pair<multimap<int,PIDs>::iterator, multimap<int,PIDs>::iterator> range;
      //   range = proc2PID.equal_range(process);
      //   for (multimap<int,PIDs>::iterator mmit=range.first; mmit!=range.second; ++mmit){
      //     int PID1 = mmit->second.PID1;
      //     int PID2 = mmit->second.PID2;
      //     double cross_section = external_xsec(PID1, PID2);       // This is the call to the external cross section tool
      //     if(cross_section == 0) cross_section = external_xsec(PID2, PID1); // Should not be necessary if external tool is sensible, my example is not
      //     NLOxsec += cross_section;
      //   }
      //   xsec.insert ( pair<int,double>(process,NLOxsec) );
      // }

        cout << DEBUG_PREFIX << "getProcessCrossSections: it = XSEC_CALCULATION, result.size() = " << result.size() << endl;          
      }

    }


    /// Compute a cross-section from Monte Carlo
    void getMCxsec(xsec& result)
    {
      using namespace Pipes::getMCxsec;

      // Don't bother if there are no analyses that will use this.
      if (Dep::RunMC->analyses.empty()) return;

      // Reset the xsec objects on all threads
      if (*Loop::iteration == COLLIDER_INIT_OMP)
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
      if (*Loop::iteration == COLLIDER_INIT_OMP)
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
      if (*Loop::iteration == COLLIDER_INIT_OMP)
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
      if (*Loop::iteration == COLLIDER_INIT_OMP)
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

    }


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
      if (*Loop::iteration == COLLIDER_INIT_OMP)
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
    }

  }
}
