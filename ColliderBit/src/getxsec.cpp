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


    /// Get cross-section info as map_str_dbl (for simple printing)
    void getXsecInfoMap(map_str_dbl& result)
    {
      using namespace Pipes::getXsecInfoMap;
      result.clear();
      const xsec& xs = (*Dep::CrossSection);
      result = xs.get_content_as_map();
    }

  }
}
