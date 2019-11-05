//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  ColliderBit module functions for calculating 
///  event weights
///
///  The weight functions in this file are
///  independent of the particular Py8Collider type
///  used in event generation.
///
///  *********************************************
///
///  Authors (add name and date if you modify):
///
///  \author Anders Kvellestad
///          (a.kvellestad@imperial.ac.uk)
///  \date   2019 Sep
///
///  *********************************************

#include "gambit/ColliderBit/ColliderBit_eventloop.hpp"

// #define COLLIDERBIT_DEBUG
#define DEBUG_PREFIX "DEBUG: OMP thread " << omp_get_thread_num() << ":  "

namespace Gambit
{

  namespace ColliderBit
  {

    extern std::map<std::string,bool> event_weight_flags;

    /// A function that sets the event weight to unity
    void _setEventWeight_unity(HEPUtils::Event& event, const BaseCollider*)  // <-- Ignoring second argument
    {
      event.set_weight(1.0);
    }

    /// Module function providing an instance of EventWeighterFunctionType
    /// pointing to _setEventWeight_unity
    void setEventWeight_unity(EventWeighterFunctionType& result)
    {
      using namespace Pipes::setEventWeight_unity;
      result = std::bind(_setEventWeight_unity, std::placeholders::_1, std::placeholders::_2);
    }



    /// A function that sets the event weight based on the process cross-sections
    void _setEventWeight_fromCrossSection(HEPUtils::Event& event, const BaseCollider* HardScatteringSim_ptr, const map_int_process_xsec& ProcessCrossSectionsMap)
    {
      // Initialize weight
      double weight = 1.0;

      // Get process code from the generator
      int process_code = HardScatteringSim_ptr->process_code();

      #ifdef COLLIDERBIT_DEBUG
        cout << DEBUG_PREFIX << "Current process_code: " << process_code << endl;
      #endif

      // Get the process_xsec_container instance that holds the externally provided cross-section for this process
      process_xsec_container xs = ProcessCrossSectionsMap.at(process_code);

      // Get the generator cross-section for this process
      double process_xsec_generator = HardScatteringSim_ptr->xsec_fb(process_code);

      #ifdef COLLIDERBIT_DEBUG
        cout << DEBUG_PREFIX << "- process_code: " << process_code << ", xsec_fb: " << HardScatteringSim_ptr->xsec_fb(process_code) << endl;
      #endif

      // Add the generator cross-sections for other process codes which also 
      // contribute to the externaly provided cross-section
      for (int other_process_code : xs.processes_sharing_xsec())
      {
        process_xsec_generator += HardScatteringSim_ptr->xsec_fb(other_process_code);
        #ifdef COLLIDERBIT_DEBUG
          cout << DEBUG_PREFIX << "- process_code: " << other_process_code << ", xsec_fb: " << HardScatteringSim_ptr->xsec_fb(other_process_code) << endl;
        #endif
      }

      // Event weight = [external cross-section] / [sum of contributing generator cross-sections]
      if (process_xsec_generator > 0.0)
      {
        weight = xs.xsec() / process_xsec_generator;
      }
      else
      {
        std::stringstream errmsg_ss;
        errmsg_ss << "Generated an event of process " << process_code << " for which the generator itself predicts a cross-section <= 0. What am I supposed to do with that?";
        ColliderBit_error().raise(LOCAL_INFO, errmsg_ss.str());
      }

      #ifdef COLLIDERBIT_DEBUG
        cout << DEBUG_PREFIX << "Total process_xsec: " << xs.xsec() << ",  process_xsec_MC: " << process_xsec_generator << ",  weight: " << weight << endl;
      #endif

      event.set_weight(weight);
    }

    /// Module function providing an instance of EventWeighterFunctionType
    /// pointing to _setEventWeight_fromCrossSection
    void setEventWeight_fromCrossSection(EventWeighterFunctionType& result)
    {
      using namespace Pipes::setEventWeight_fromCrossSection;

      static bool first = true;
      if (first)
      {
        event_weight_flags["weight_by_cross_section"] = true;
        first = false;
      }

      if(*Loop::iteration < 0) return;

      result = std::bind(_setEventWeight_fromCrossSection,
                         std::placeholders::_1,
                         std::placeholders::_2,
                         *Dep::ProcessCrossSectionsMap);
    }


  } 
} 


