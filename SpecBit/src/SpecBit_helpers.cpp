//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  Definitions of convenience (i.e. non-Gambit)
///  functions used by more than one SpecBit 
///  source file.
///
///  *********************************************
///
///  Authors (add name and date if you modify):
///
///  \author Tomas Gonzalo
///          (tomas.gonzalo@monash.edu)
///  \date 2019 Nov
///
///  *********************************************

#include "gambit/SpecBit/SpecBit_helpers.hpp"

namespace Gambit
{

  namespace SpecBit
  {

    /// @{

    /// Non-Gambit helper functions
    //  =======================================================================
    //  These are not known to Gambit, but perform helper tasks used by the
    //  Gambit module functions.

    double get_b_pole(SMInputs sminputs)
    {
      // TODO: Should do some computations here to return the pole mass
      return sminputs.mBmB;
    }

  }
}

