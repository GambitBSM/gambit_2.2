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

      // reset all bounceActionThreshold (non-thermal and thermal) entries to -1
      static std::vector<std::string> bounceAction = {"bounceActionThreshold", "bounceActionThresholdThermal"};
      // maximum length of bounce action vector with currenct vevacious interface
      static int max_length = 4;

      // loop through thermal/non-thermal entries & vector length
      for(auto&& action_str: bounceAction)
      {
        for(int ii =0; ii < max_length; ii ++)
        {
          result_map[panic_vacuum][ action_str+"_["+std::to_string(ii)+"]"] = -1;  
        }
      }
      
    }

    /// add a SpectrumEntry type to the 'spec_entry_map' map. GAMBIT will iterate through it and 
    /// pass all contents to vevacious before it is called.
    void SpectrumEntriesForVevacious::add_entry (std::string name, vec_pair_int_dbl parameters, int dimension)
    { 
        // create spectrum entry 
        SpectrumEntry entry;
        entry.name = name;
        entry.parameters = parameters;
        entry.dimension = dimension;

        // .. and store it in map
        spec_entry_map[entry.name]= entry;
    };
  }
}
