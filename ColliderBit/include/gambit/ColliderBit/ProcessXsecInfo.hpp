//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  ColliderBit class for holding the cross-section
///  (xsec instance) for a given process code, together
///  info needed for correct event weighting
///
///  *********************************************
///
///  Authors (add name and date if you modify):
///
///  \author Anders Kvellestad
///          (a.kvellestad@imperial.ac.uk)
///  \date 2019 Sept
///
///  *********************************************

#include <vector>

#include "gambit/ColliderBit/xsec.hpp"

#pragma once

namespace Gambit
{

  namespace ColliderBit
  {

    /// A struct for holding a xsec instance and associated info 
    /// for a process identified by the collider process code
    struct ProcessXsecInfo
    {
        typedef std::pair<int,int> PID_pair;
        typedef std::vector<std::pair<int,int>> vec_PID_pairs;

        int process_code;
        vec_PID_pairs pid_pairs;
        xsec process_xsec;
        std::vector<int> processes_sharing_xsec;
    };

  }
}
