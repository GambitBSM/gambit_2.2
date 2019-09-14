//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  ColliderBit module functions for calculating 
///  event weights
/// 
///  The weight functions in this file are
///  specific to the oarticular Py8Collider type
///  used for event generation.
///
///  *********************************************
///
///  Authors (add name and date if you modify):
///
///  \author Anders Kvellestad
///          (a.kvellestad@imperial.ac.uk)
///  \date   2019 Sept
///
///  *********************************************

#include "gambit/ColliderBit/ColliderBit_eventloop.hpp"

#define COLLIDERBIT_DEBUG
#define DEBUG_PREFIX "DEBUG: OMP thread " << omp_get_thread_num() << ":  "

namespace Gambit
{

  namespace ColliderBit
  {

    /// A function that sets the event weight based on the process cross-sections
    template<typename PythiaT, typename EventT>
    void _setEventWeight_fromCrossSection(HEPUtils::Event& event, const BaseCollider* collider_ptr, const map_int_ProcessXsecInfo& ProcessCrossSections)
    {
      cout << DEBUG_PREFIX << ": _setEventWeight_fromCrossSection: Starting function..." << endl;

      double weight = 1.0;

      // Cast the collider ptr to the correct type
      const Py8Collider<PythiaT,EventT>* HardScatteringSim_ptr = dynamic_cast<const Py8Collider<PythiaT,EventT>*>(collider_ptr);

      // Get process code from Pythia
      int process_code = HardScatteringSim_ptr->pythia()->info.code();

      // Get the ProcessXsecInfo instance that holds the externally provided cross-section for this process
      ProcessXsecInfo xs_info = ProcessCrossSections.at(process_code);

      // Pythia cross-section for this process
      double process_xsec_pythia = HardScatteringSim_ptr->pythia()->info.sigmaGen(process_code) * 1e-9;  // Pythia uses mb, we use pb

      cout << std::scientific << std::setprecision(5);
      cout << DEBUG_PREFIX << ": info.sigmaGen(" << process_code << "): " << HardScatteringSim_ptr->pythia()->info.sigmaGen(process_code) * 1e-9 << endl;

      // Add the Pythia cross-sections for other process codes which also 
      // contribute to the externaly provided cross-section
      for (int other_process_code : xs_info.processes_sharing_xsec)
      {
        process_xsec_pythia += HardScatteringSim_ptr->pythia()->info.sigmaGen(other_process_code) * 1e-9;  // Pythia uses mb, we use pb
        cout << DEBUG_PREFIX << ": info.sigmaGen(" << other_process_code << "): " << HardScatteringSim_ptr->pythia()->info.sigmaGen(other_process_code) * 1e-9 << endl;
      }

      // Event weight = external cross-section / sum of contributing Pythia cross-sections
      if (process_xsec_pythia > 0.0)
      {
        weight = xs_info.process_xsec() / process_xsec_pythia;
      }
      else
      {
        std::stringstream errmsg_ss;
        errmsg_ss << "Pythia generated an event of process " << process_code << " for which itself predicts a cross-section <= 0.0 pb... What am I supposed to do with that?";
        ColliderBit_error().raise(LOCAL_INFO, errmsg_ss.str());
      }

      cout << std::scientific << std::setprecision(5);
      cout << DEBUG_PREFIX << "process_code: " << process_code << ",  process_xsec: " << xs_info.process_xsec() << ",  process_xsec_pythia: " << process_xsec_pythia << ",  weight: " << weight << endl;

      event.set_weight(weight);
    }


    /// Construct a module function that provides an instance of 
    /// EventWeighterType_Py8Collider pointing to _setEventWeight_fromCrossSection
    #define SET_EVENT_WEIGHT_FROM_CROSS_SECTION(NAME, PYTHIA_NS)                    \
    void NAME(EventWeighterType_Py8Collider& result)                                \
    {                                                                               \
      using namespace Pipes::NAME;                                                  \
                                                                                    \
      if(*Loop::iteration < 0) return;                                              \
                                                                                    \
      result = std::bind(_setEventWeight_fromCrossSection<PYTHIA_NS::Pythia8::Pythia, PYTHIA_NS::Pythia8::Event>,  \
                         std::placeholders::_1,                                     \
                         std::placeholders::_2,                                     \
                         *Dep::ProcessCrossSections);                               \
    }

  } 
} 


