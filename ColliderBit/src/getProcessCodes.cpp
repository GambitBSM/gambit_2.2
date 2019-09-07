//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  ColliderBit module functions for dealing with
///  lists of the active processes
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


    // enum specialIterations { BASE_INIT = -1,
    //                          COLLIDER_INIT = -2,
    //                          START_SUBPROCESS = -3,
    //                          COLLECT_CONVERGENCE_DATA = -4,
    //                          CHECK_CONVERGENCE = -5,
    //                          END_SUBPROCESS = -6,
    //                          COLLIDER_FINALIZE = -7,
    //                          BASE_FINALIZE = -8};


    /// Get the list of active Pythia process codes
    void getPythiaProcessCodes(std::vector<int>& result)
    {
      using namespace Pipes::getPythiaProcessCodes;

      if (*Loop::iteration == COLLIDER_INIT)
      {
        result.clear();
      }

      if (*Loop::iteration == XSEC_CALCULATION)
      {
        result = Dep::HardScatteringSim->codesHard();
        // result = std::vector<int>({1,2,3,4,5});
      }

      // if (*Loop::iteration == BASE_INIT) { cout << DEBUG_PREFIX << "getPythiaProcessCodes: in BASE_INIT: result.size() = " << result.size() << endl; }
      // if (*Loop::iteration == COLLIDER_INIT) { cout << DEBUG_PREFIX << "getPythiaProcessCodes: in COLLIDER_INIT: result.size() = " << result.size() << endl; }
      // if (*Loop::iteration == COLLIDER_INIT_OMP) { cout << DEBUG_PREFIX << "getPythiaProcessCodes: in COLLIDER_INIT_OMP: result.size() = " << result.size() << endl; }
      // if (*Loop::iteration == XSEC_CALCULATION) { cout << DEBUG_PREFIX << "getPythiaProcessCodes: in XSEC_CALCULATION: result.size() = " << result.size() << endl; }
      // if (*Loop::iteration == START_SUBPROCESS) { cout << DEBUG_PREFIX << "getPythiaProcessCodes: in START_SUBPROCESS: result.size() = " << result.size() << endl; }
      // if (*Loop::iteration == 0) { cout << DEBUG_PREFIX << "getPythiaProcessCodes: in ITERATION 0: result.size() = " << result.size() << endl; }
      // if (*Loop::iteration == 1) { cout << DEBUG_PREFIX << "getPythiaProcessCodes: in ITERATION 1: result.size() = " << result.size() << endl; }
      // if (*Loop::iteration == END_SUBPROCESS) { cout << DEBUG_PREFIX << "getPythiaProcessCodes: in END_SUBPROCESS: result.size() = " << result.size() << endl; }
      // if (*Loop::iteration == COLLIDER_FINALIZE) { cout << DEBUG_PREFIX << "getPythiaProcessCodes: in COLLIDER_FINALIZE: result.size() = " << result.size() << endl; }
      // if (*Loop::iteration == BASE_FINALIZE) { cout << DEBUG_PREFIX << "getPythiaProcessCodes: in BASE_FINALIZE: result.size() = " << result.size() << endl; }

    }


    /// Translate list of Pythia process codes to list of (PID,PID) pairs
    /// for the two final state particles of the hard process

    // void getProcessPIDPairs(std::vector<std::pair<int,int>>& result)
    // void getProcessPIDPairs(std::vector<PID_pair>& result)
    void getProcessPIDPairs(vec_PID_pairs& result)
    {
      using namespace Pipes::getProcessPIDPairs;

      if (*Loop::iteration == COLLIDER_INIT)
      {
        result.clear();
      }

      if (*Loop::iteration == XSEC_CALCULATION)
      {
        std::vector<int> process_codes = *Dep::ProcessCodes;
        for (int c : process_codes)
        {
          result.push_back( PID_pair(c,c) );
        }
      }

      // _Anders
      if (*Loop::iteration == START_SUBPROCESS)
      {
        cout << DEBUG_PREFIX << "getProcessPIDPairs: it = START_SUBPROCESS, result.size() = " << result.size() << endl;
      }

    }

  } 
} 


