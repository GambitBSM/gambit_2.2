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
        // constructor initialises every member to -1 to avoid 
        // problems when printing results when vevacious did not run
        VevaciousResultContainer();
        ~VevaciousResultContainer();

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
    };
  }
}

#endif // defined __SpecBit_types_hpp__