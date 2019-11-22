//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  Declarations of convenience (i.e. non-Gambit)
///  functions used by more than one SpecBit 
///  source file.
///
///  *********************************************
///
///  Authors (add name and date if you modify):
///
///  \author Ben Farmer
///          (benjamin.farmer@fysik.su.se)
///    \date 2014 Sep - Dec, 2015 Jan - May
///  
///  \author Tomas Gonzalo
///          (tomas.gonzalo@monash.edu)
///  \date 2019 Nov
///
///  *********************************************

#ifndef __SpecBit_helpers_hpp__
#define __SpecBit_helpers_hpp__

#include "gambit/Elements/gambit_module_headers.hpp"
#include "gambit/Elements/sminputs.hpp"
#include "gambit/Elements/spectrum.hpp"

#include "gambit/SpecBit/SpecBit_rollcall.hpp"

#include "gambit/Models/partmap.hpp"

namespace Gambit
{

  namespace SpecBit
  {

    /// @{

    /// Non-Gambit helper functions
    //  =======================================================================
    //  These are not known to Gambit, but perform helper tasks used by the
    //  Gambit module functions.

    /// Check that the spectrum has the canonical LSP for the model being scanned.
    void check_LSP(const Spectrum& spec, std::vector<int> LSPs);

    /// Helper to work with pointer
    void check_LSP(const Spectrum* spec, std::vector<int> LSPs);

    /// Add gravitino mass to the spectrum and list of LSPs
    void add_gravitino_mass(Spectrum& spec, std::vector<int> &LSPs, double mG, const safe_ptr<Options>& runOptions);
  



  }
}
 
#endif
