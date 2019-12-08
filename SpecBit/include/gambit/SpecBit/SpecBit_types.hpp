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
#include "gambit/Utils/util_types.hpp"

namespace Gambit
{

  namespace SpecBit
  {

    // thisis the type we need to pass the spectrum entries to vevacious
    typedef std::vector<std::pair<int,double>> vec_pair_int_dbl;

    // each of the vectors defined above have another int-type setting that 
    // needs to be passed to vevacious to read in the vector (don't know what it is...)
    // so we need to safe this information as well to be able to read in the stuff properly
    typedef std::pair<int,vec_pair_int_dbl> pair_int_spectrum_vec;

    // now the map to this horrible type
    typedef std::map<std::string,pair_int_spectrum_vec> map_str_pair_int_spectrum_vec;


    // => in the end we want a map looking something like this
    // map[name] = [int other_setting, vector_pair_int_dbl]
    class SpectrumEntriesForVevacious
    {
        public:
            SpectrumEntriesForVevacious(){};

            // setter functions for scale, inputPath and input Filenames
            void set_scale (double inScale) {scale = inScale;};
            void set_inputPath (std::string inPath) {inputPath = inPath;};
            void set_inputFilename (std::string inFile) {inputFilename = inFile;};
            
            // getter functions for scale, inputPath and input Filenames
            double get_scale () {return scale;};
            std::string get_inputFilename () {return inputFilename;};
            std::string get_inputPath () {return inputPath;};

            // adds an entry to the spec_entry_map
            void add_entry (std::string name, vec_pair_int_dbl vec, int setting);

            // return spec_entry_map -> iterate through it to pass all entries to vevacious
            map_str_pair_int_spectrum_vec get_spec_entry_map() {return spec_entry_map;};

        private:
            double scale;
            int other_setting;
            std::string name;
            std::string inputFilename;
            std::string inputPath;
            vec_pair_int_dbl vec;
            map_str_pair_int_spectrum_vec spec_entry_map;

    };


    /* Class that stores the results computed by vevacious that will be 
    needed by other capabilites in GAMBIT  */
    class VevaciousResultContainer
    {

      public:
        // constructor initialises every member to -1 to avoid 
        // problems when printing results when vevacious did not run
        VevaciousResultContainer();

        // reset all memember variables to value -1
        void reset_results();
        void vevacious_ran(){vevaciousRunFlag = true;};
        void vevacious_ran_reset(){vevaciousRunFlag = false;};

        // setter functions for members, set bool thermal to true
        // to set thermal values
        void set_lifetime(double val) {lifetime = val;};
        void set_thermalProbability(double val) {thermalProbability = val;};
        
        void set_bounceActionThreshold  (double val, bool thermal);
        void set_bounceActionStraight   (double val, bool thermal);
        void set_firstPathFinder        (double val, bool thermal);
        void set_secondPathFinder       (double val, bool thermal);

        // getter functions for members, set thermal to ture to 
        // get thermal values
        double get_lifetime()           {return lifetime;};
        double get_thermalProbability() {return thermalProbability;};

        double get_bounceActionThreshold(bool thermal);
        double get_bounceActionStraight (bool thermal);
        double get_firstPathFinder      (bool thermal);
        double get_secondPathFinder     (bool thermal);

        void set_results (str panic_vaccum, str name, double val){result_map[panic_vaccum][name]=val;}
        void clear_results(str panic_vaccum);

        map_str_dbl get_nearest_results() {return result_map["nearest"];}
        map_str_dbl get_global_results() {return result_map["global"];}
        double get_lifetime(str panic_vaccum) {return result_map[panic_vaccum]["lifetime"];};
        
      private:
        bool vevaciousRunFlag; 
        double lifetime;
        double thermalProbability;

        double bounceActionThreshold;
        double bounceActionThresholdThermal;

        double bounceActionStraight;
        double bounceActionStraightThermal;

        double firstPathFinder;
        double firstPathFinderThermal;

        double secondPathFinder;
        double secondPathFinderThermal;  

        //map_str_dbl global_results;
        //map_str_dbl nearest_results;   
        map_str_map_str_dbl result_map;   
    };
  }
}

#endif // defined __SpecBit_types_hpp__