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
///  \author Anders Kvellestad
///          (a.kvellestad@imperial.ac.uk)
///  \date 2019 Sep
///
///  *********************************************

#include "gambit/ColliderBit/ColliderBit_eventloop.hpp"

#define COLLIDERBIT_DEBUG
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

      return std::make_pair(xsec_fb, xsec_uncert_fb);
    }



    // ======= Module functions =======


    /// Dummy function for testing out how to return a PIDPairCrossSectionFunc
    xsec PIDPairCrossSection_dummy(const PID_pair& pids, const Spectrum& MSSM_spectrum)
    {
      xsec xs_result;

      // Get an SLHA1 object
      const SLHAstruct& slha = MSSM_spectrum.getSLHAea(1);

      // "Calculate" cross section
      // (only stop-stopbar production)
      double xs_fb = 0.0;
      if (pids.first == 1000002 && pids.second == -1000002)
      {
        xs_fb = 3.09816e-11 + 9.08223e-13;
      }
      else
      {
        xs_fb = 1.0;
      }
      double xs_rel_err = 0.01;
      // double xs_err = xs_fb * xs_rel_err;

      // Save result in xs_result
      xs_result.set_xsec(xs_fb, xs_rel_err);

      // Construct info string of the form "PID1:<PID1>, PID2:<PID2>"
      std::stringstream info_ss;
      info_ss << "PID1:" << pids.first << ", " << "PID2:" << pids.second;
      xs_result.set_info_string(info_ss.str());

      return xs_result;
    }

    /// Get a PIDPairCrossSectionFunc pointing to PIDPairCrossSection_dummy,
    /// with the second argument already filled by *Dep::MSSM_spectrum
    void getPIDPairCrossSectionFunc_dummy(PIDPairCrossSectionFuncType& result)
    {
      using namespace Pipes::getPIDPairCrossSectionFunc_dummy;
      result = std::bind(PIDPairCrossSection_dummy, std::placeholders::_1, *Dep::MSSM_spectrum);
    }



    /// Get the cross-section from the xsec_example backend
    xsec PIDPairCrossSection_xsec_example(PID_pair pids, 
                                          const Spectrum& MSSM_spectrum)
                                          // double (*xsec_fb_fptr)(PID_pair&, map_str_dbl&, map_str_bool&))
    {
      xsec xs_result;

      // Get an SLHA1 object
      const SLHAstruct& slha = MSSM_spectrum.getSLHAea(1);

      // Call cross-section calculation 
      map_str_dbl proc_params;
      map_str_bool proc_flags;

      double xs_fb = 3.1415;

      // double xs_fb = 1.0e-3 * xsec_fb_fptr(pids, proc_params, proc_flags);

      double xs_rel_err = 0.01;

      // Save result in xs_result
      xs_result.set_xsec(xs_fb, xs_rel_err);

      // Construct info string of the form "PID1:<PID1>, PID2:<PID2>"
      std::stringstream info_ss;
      info_ss << "PID1:" << pids.first << ", " << "PID2:" << pids.second;
      xs_result.set_info_string(info_ss.str());

      return xs_result;
    }


    /// Get a PIDPairCrossSectionFunc pointing to the xsec_example backend
    void getPIDPairCrossSectionFunc_xsec_example(PIDPairCrossSectionFuncType& result)
    {
      using namespace Pipes::getPIDPairCrossSectionFunc_xsec_example;

      // result = std::bind(PIDPairCrossSection_xsec_example, std::placeholders::_1, *Dep::MSSM_spectrum, BEreq::xsec_example_xsec_fb.pointer());
      result = std::bind(PIDPairCrossSection_xsec_example, std::placeholders::_1, *Dep::MSSM_spectrum);

      // Bind to the backend function xsec_fb with the two first arguments fixed

      // // Get an SLHA1 object
      // const SLHAstruct& slha = MSSM_spectrum.getSLHAea(1);

      // _Anders
      // BACKEND_REQ(xsec_example_xsec_fb, (), double, (std::vector<double>&, map_str_dbl&, map_str_bool&))

      // map_str_dbl proc_params;
      // map_str_bool proc_flags;

      // result = std::bind(BEreq::xsec_example_xsec_fb.pointer(), std::placeholders::_1, proc_params, proc_flags);
    }






    /// Get a map between Pythia process codes and cross-sections
    void getProcessCrossSectionsMap(map_int_ProcessXsecInfo& result)
    {
      using namespace Pipes::getProcessCrossSectionsMap;

      // Use a static variable to communicate the result calculated on thread 0 during 
      // iteration XSEC_CALCULATION to all threads during iteration START_SUBPROCESS
      static map_int_ProcessXsecInfo shared_result;

      // Only thread 0
      if(*Loop::iteration == COLLIDER_INIT)
      {
        shared_result.clear();
      }

      // All threads
      if(*Loop::iteration == COLLIDER_INIT_OMP)
      {
        result.clear();
      }

      // Only thread 0
      if(*Loop::iteration == XSEC_CALCULATION)
      {
        cout << DEBUG_PREFIX << "getProcessCrossSectionsMap: it = XSEC_CALCULATION, ProcessCodes.size() = " << Dep::ProcessCodes->size() << endl;          
        cout << DEBUG_PREFIX << "getProcessCrossSectionsMap: it = XSEC_CALCULATION, ProcessCodeToPIDPairsMap.size() = " << Dep::ProcessCodeToPIDPairsMap->size() << endl;          

        // Loop over all active processes and construct the cross-section map (shared_result)
        for (size_t i = 0; i != Dep::ProcessCodes->size(); ++i)
        {
          // Get process code
          int current_pcode = Dep::ProcessCodes->at(i);

          // Construct a ProcessXsecInfo instance to be stored in the shared_result map
          ProcessXsecInfo xs_info;
          xs_info.process_code = current_pcode;

          // Get iterator bounds (as a pair) over the multimap entries that match the key current_pcode
          auto mm_range = Dep::ProcessCodeToPIDPairsMap->equal_range(current_pcode);

          // Loop over these elements in the multimap
          for (auto mm_it = mm_range.first; mm_it != mm_range.second; ++mm_it)
          {
            const PID_pair& pids = mm_it->second;

            // Call cross-section calculator here.
            // (We're gonna assume that the calculator is smart enough to re-order two PIDs if it needs to.)
            xsec xs = (*Dep::PIDPairCrossSectionFunc)(pids);

            // Accumulate result in the ProcessXsecInfo::process_xsec variable
            xs_info.process_xsec.sum_xsecs(xs);
            xs_info.pid_pairs.push_back(pids);
          }

          // Construct info string of the form "ProcessCode:<current_pcode>"
          std::stringstream info_ss;
          info_ss << "ProcessCode:" << current_pcode;
          xs_info.process_xsec.set_info_string(info_ss.str());


          // Now we figure out if the current_pcode process shares the cross-section
          // stored in in xs_info.process_xsec with any other process codes

          // Loop over *all* entries (process code <--> PID pair) in the multimap Dep::ProcessCodeToPIDPairsMap
          for (auto mm_it = Dep::ProcessCodeToPIDPairsMap->begin(); mm_it != Dep::ProcessCodeToPIDPairsMap->end(); ++mm_it)
          {
            // Extract the process code (pc) and PID pair (pp)
            int pc = mm_it->first;
            const PID_pair& pp = mm_it->second;

            if (pc == current_pcode) continue;

            // @todo What's the right choice here?
            // // Don't add more copies of the same process code! ...Or should we?
            // if(std::find(xs_info.processes_sharing_xsec.begin(), xs_info.processes_sharing_xsec.end(), pc) != xs_info.processes_sharing_xsec.end()) 

            // Check if the PID pair pp mathces one of the PID pairs for the current_pcode process
            if(std::find(xs_info.pid_pairs.begin(), xs_info.pid_pairs.end(), pp) != xs_info.pid_pairs.end()) 
            {
              // Check that pc is itself in one of the active processes, i.e. listed in Dep::ProcessCodes
              if(std::find(Dep::ProcessCodes->begin(), Dep::ProcessCodes->end(), pc) != Dep::ProcessCodes->end())  
              {
                // Add pc to the list of processes that share cross-section with current_pcode
                xs_info.processes_sharing_xsec.push_back(pc);
              }
              else
              {
                std::stringstream errmsg_ss;
                errmsg_ss << "For correct cross-section scaling of collider process " << current_pcode;
                errmsg_ss << ", process " << pc << " must also be activated. Please check your collider settings.";
                ColliderBit_error().raise(LOCAL_INFO, errmsg_ss.str());
              }
            }
          }

          // Store xs_info in the shared_result map
          shared_result[current_pcode] = xs_info;
        }

        // Let thread 0 return the correct result already after iteration XSEC_CALCULATION
        result = shared_result;
      }


      // All threads
      if (*Loop::iteration == START_SUBPROCESS)
      {
        // All threads read the result from shared_result
        result = shared_result;
      }

    }


    /// Compute a cross-section from Monte Carlo
    void getMCxsec(MC_xsec& result)
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
          const double xs_fb = (*Dep::HardScatteringSim)->xsec_fb();
          const double xserr_fb = (*Dep::HardScatteringSim)->xsecErr_fb();
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

    /// Return MC_xsec as const base_xsec*
    void get_MC_xsec_as_base(const base_xsec*& result)
    {
      using namespace Pipes::get_MC_xsec_as_base;
      result = &(*Dep::TotalCrossSectionFromMC);
    }

    /// Return xsec as const base_xsec*
    void get_xsec_as_base(const base_xsec*& result)
    {
      using namespace Pipes::get_xsec_as_base;
      result = &(*Dep::TotalCrossSection);
    }



    /// Get a cross-section from NLL-FAST
    void getNLLFastxsec(xsec& result)
    {
      using namespace Pipes::getNLLFastxsec;

      // Use a static variable to communicate the result calculated on thread 0 during 
      // iteration XSEC_CALCULATION to all threads during iteration START_SUBPROCESS
      static xsec shared_result;

      // Don't bother if there are no analyses that will use this.
      if (Dep::RunMC->analyses.empty()) return;

      // Only thread 0
      if(*Loop::iteration == COLLIDER_INIT) shared_result.reset();
      
      // All threads
      if (*Loop::iteration == COLLIDER_INIT_OMP) result.reset();

      // Set the xsec and its error.
      // Only thread 0
      if (*Loop::iteration == XSEC_CALCULATION)
      {
        double xs_fb = 0.1;             // replace with xsec from NLL-Fast
        double xserr_fb = 0.1 * xs_fb;  // or whatever
        shared_result.set_xsec(xs_fb, xserr_fb);

        // Let thread 0 return the correct result already after iteration XSEC_CALCULATION
        result = shared_result;
      }

      // All threads
      if (*Loop::iteration == START_SUBPROCESS)
      {
        // All threads copy the result from shared_result
        result = shared_result;
      }

    }


    /// A function that reads the total cross-section from the input file, but builds up the number of events from the event loop
    void getYAMLxsec(xsec& result)
    {
      using namespace Pipes::getYAMLxsec;

      // Use a static variable to communicate the result calculated on thread 0 during 
      // iteration XSEC_CALCULATION to all threads during iteration START_SUBPROCESS
      static xsec shared_result;

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

      // Only thread 0
      if(*Loop::iteration == COLLIDER_INIT) shared_result.reset();
      
      // All threads
      if (*Loop::iteration == COLLIDER_INIT_OMP) result.reset();

      // Set the xsec and its error
      // Only thread 0
      if (*Loop::iteration == XSEC_CALCULATION)
      {
        std::pair<double,double> temp = convert_xsecs_to_fb(input_xsec, input_xsec_uncert, input_unit, input_fractional_uncert);
        double xsec_fb = temp.first;
        double xsec_uncert_fb = temp.second;
        shared_result.set_xsec(xsec_fb, xsec_uncert_fb);

        // Let thread 0 return the correct result already after iteration XSEC_CALCULATION
        result = shared_result;
      }

      // All threads
      if (*Loop::iteration == START_SUBPROCESS)
      {
        // All threads copy the result from shared_result
        result = shared_result;
      }

    }



    /// A function that reads a list of (SLHA file, total cross-section) pairs from the input YAML file
    void getYAMLxsec_SLHA(xsec& result)
    {
      using namespace Pipes::getYAMLxsec_SLHA;

      // Use a static variable to communicate the result calculated on thread 0 during 
      // iteration XSEC_CALCULATION to all threads during iteration START_SUBPROCESS
      static xsec shared_result;

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

      // Only thread 0
      if(*Loop::iteration == COLLIDER_INIT) shared_result.reset();
      
      // All threads
      if (*Loop::iteration == COLLIDER_INIT_OMP) result.reset();

      // Set the xsec and its error
      if (*Loop::iteration == XSEC_CALCULATION)
      {
        double input_xsec = colOptions_xsec.getValue<double>(filename);
        double input_xsec_uncert = colOptions_uncert.getValue<double>(filename);

        std::pair<double,double> temp = convert_xsecs_to_fb(input_xsec, input_xsec_uncert, input_unit, input_fractional_uncert);
        double xsec_fb = temp.first;
        double xsec_uncert_fb = temp.second;
        shared_result.set_xsec(xsec_fb, xsec_uncert_fb);

        // Let thread 0 return the correct result already after iteration XSEC_CALCULATION
        result = shared_result;
      }

      // All threads
      if (*Loop::iteration == START_SUBPROCESS)
      {
        // All threads copy the result from shared_result
        result = shared_result;
      }

    }


    /// A function that assigns a total cross-sections directly from the scan parameters
    /// (for models CB_SLHA_simpmod_scan_model and CB_SLHA_scan_model)
    void getYAMLxsec_param(xsec& result)
    {
      using namespace Pipes::getYAMLxsec_param;

      // Use a static variable to communicate the result calculated on thread 0 during 
      // iteration XSEC_CALCULATION to all threads during iteration START_SUBPROCESS
      static xsec shared_result;

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

      // Only thread 0
      if(*Loop::iteration == COLLIDER_INIT) shared_result.reset();
      
      // All threads
      if (*Loop::iteration == COLLIDER_INIT_OMP) result.reset();

      // Set the xsec and its error
      // Only thread 0
      if (*Loop::iteration == XSEC_CALCULATION)
      {
        double input_xsec = *Param.at(xsec_pnames.first);
        double input_xsec_uncert = *Param.at(xsec_pnames.second); 

        std::pair<double,double> temp = convert_xsecs_to_fb(input_xsec, input_xsec_uncert, input_unit, input_fractional_uncert);
        double xsec_fb = temp.first;
        double xsec_uncert_fb = temp.second;
        shared_result.set_xsec(xsec_fb, xsec_uncert_fb);

        // Let thread 0 return the correct result already after iteration XSEC_CALCULATION
        result = shared_result;
      }

      // All threads
      if (*Loop::iteration == START_SUBPROCESS)
      {
        // All threads copy the result from shared_result
        result = shared_result;
      }

    }


    /// Get cross-section info as map_str_dbl (for simple printing)
    void getTotalCrossSectionAsMap(map_str_dbl& result)
    {
      using namespace Pipes::getTotalCrossSectionAsMap;

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
        for(auto s_d_pair : (*Dep::TotalCrossSection)->get_content_as_map())
        {
          std::string new_key(Dep::RunMC->current_collider());
          new_key.append("__").append(s_d_pair.first);
          result[new_key] = s_d_pair.second;
        }
      }
    }

  }
}
