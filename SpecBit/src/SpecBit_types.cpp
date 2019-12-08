//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  Source code for types for module SpecBit.
///  For instructions on adding new types, see
///  the corresponding header.
///
///  *********************************************
///
///  Authors (add name and date if you modify):
///
///  \author Janina Renk
///          (janina.renk@fysik.su.se)
///  \date 2019 July
///  *********************************************

#include <string>
#include <iostream>
#include <stdlib.h>     /* malloc, free, rand */
#include <valarray>

#include "gambit/SpecBit/SpecBit_types.hpp"

namespace Gambit
{
  namespace SpecBit
  {

    /// delete all entries of the vevacious results map and set
    /// lifetime, probability and bounceActionThreshould to -1
    /// to be able to easily filter point out for which vevacious
    /// did not run successfully 
    void VevaciousResultContainer::clear_results(str panic_vacuum)
    { 
      result_map[panic_vacuum].clear();
      std::cout << " VevaciousResultContainer clearing all entries for " << panic_vacuum << std::endl;
      result_map[panic_vacuum]["lifetime"] = -1;                
      result_map[panic_vacuum]["thermalProbability"] = -1;
      result_map[panic_vacuum]["bounceActionThreshold_[0]"] = -1;  
      result_map[panic_vacuum]["bounceActionThresholdThermal_[0]"] = -1;

    }

    void SpectrumEntriesForVevacious::add_entry (std::string name, vec_pair_int_dbl vec, int setting)
    {
        pair_int_spectrum_vec pair = make_pair(setting, vec);
        spec_entry_map[name]= pair;
    };
  }
}
