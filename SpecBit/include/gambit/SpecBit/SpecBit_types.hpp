//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  Type definition header for module SpecBit.
///
///  Compile-time registration of type definitions
///  required for the rest of the code to
///  communicate with SpecBit.
///
///  Add to this if you want to define a new type
///  for the functions in SpecBit to return, but
///  you don't expect that type to be needed by
///  any other modules.
///
///  *********************************************
///
///  Authors (add name and date if you modify):
///
///  \author Janina Renk
///          (janina.renk@fysik.su.se)
///  \date 2019 July
///
///  *********************************************


#ifndef __SpecBit_types_hpp__
#define __SpecBit_types_hpp__

#include <string>
#include <iostream>
#include <vector>
#include <stdlib.h>     /* malloc, free, rand */

namespace Gambit
{

  namespace SpecBit
  {


    /* Class that stores the results computed by vevacious that will be 
    needed by other capabilites in GAMBIT  */
    class VevaciousResultContainer
    {

      public:
        VevaciousResultContainer(){PathFinderResults.clear();PathFinderThermalResults.clear();};
        ~VevaciousResultContainer(){std::cout << "AAHHHGHG I DIED" << std::endl; };

        // Methods to add entry to vevacious_result_map
        void addEntry(str key, double value) {vevacious_result_map[key]=value;};
        // return 1 if key found, zero otherwise
        int hasKey(str key){return vevacious_result_map.count(key);};
        // clears all entries from vevacious_result_map
        void clear(){vevacious_result_map.clear();};
        

        // add Entries to the PathFinderResults 
        void addPathFinderEntry         (double value) {PathFinderResults.push_back(value);};
        void addPathFinderThermalEntry  (double value) {PathFinderThermalResults.push_back(value);};

        // return PathFinderResult vectors
        std::vector<double> return_PathFinderResults(){return PathFinderResults;};
        std::vector<double> return_PathFinderThermalResults(){return PathFinderThermalResults;};
        map_str_dbl return_result_map(){return vevacious_result_map;};


      private:
        map_str_dbl vevacious_result_map;
        std::vector<double> PathFinderResults;
        std::vector<double> PathFinderThermalResults;
        
    };
  }
}

#endif // defined __SpecBit_types_hpp__