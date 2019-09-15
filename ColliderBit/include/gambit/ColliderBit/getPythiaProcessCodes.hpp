//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  ColliderBit module function for obtaining the 
///  list of the active collider processes
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

    /// Extract the list of active Pythia process codes from a Py8Collider type
    template<typename PythiaT, typename EventT>
    void _getPythiaProcessCodes(std::vector<int>& process_codes,
                                const Py8Collider<PythiaT,EventT>& HardScatteringSim,
                                const int iteration)
    {
      if (iteration == COLLIDER_INIT)
      {
        process_codes.clear();
      }

      if (iteration == XSEC_CALCULATION)
      {
        process_codes = HardScatteringSim.codesHard();
      }
    }


    /// Construct a module function providing a list of active processes
    /// for a specific Pythia (a specific Py8Collider specialization)
    #define GET_PYTHIA_PROCESS_CODES(NAME)                       \
    void NAME(std::vector<int>& result)                          \
    {                                                            \
      using namespace Pipes::NAME;                               \
      _getPythiaProcessCodes(result,                             \
                             *Dep::HardScatteringSim,            \
                             *Loop::iteration);                  \
    }



    // /// Get the list of active Pythia process codes
    // void getPythiaProcessCodes(std::vector<int>& result)
    // {
    //   using namespace Pipes::getPythiaProcessCodes;

    //   if (*Loop::iteration == COLLIDER_INIT)
    //   {
    //     result.clear();
    //   }

    //   if (*Loop::iteration == XSEC_CALCULATION)
    //   {
    //     result = Dep::HardScatteringSim->codesHard();
    //   }
    // }


  } 
} 


