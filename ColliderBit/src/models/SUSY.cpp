//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  SUSY-specific sources for ColliderBit.
///
///  *********************************************
///
///  Authors (add name and date if you modify):
///
///  \author Pat Scott
///          (p.scott@imperial.ac.uk)
///  \date 2019 Jan
///
///  *********************************************

#include "gambit/ColliderBit/getPy8Collider.hpp"
#include "gambit/ColliderBit/generateEventPy8Collider.hpp"

namespace Gambit
{
  namespace ColliderBit
  {

    // Get Monte Carlo event generator
    GET_SPECIFIC_PYTHIA(getPythia, Pythia_default, MSSM_spectrum, , IS_SUSY)
    GET_PYTHIA_AS_BASE_COLLIDER(getPythiaAsBase)

    // Run event generator
    GET_PYTHIA_EVENT(generateEventPythia)


    // Get Monte Carlo event generator from SLHA file input
    GET_SPECIFIC_PYTHIA_SLHA(getPythia_SLHA, Pythia_default, )
    // GET_PYTHIA_AS_BASE_COLLIDER(getPythia_SLHAAsBase)

    // Run event generator
    GET_PYTHIA_EVENT(generateEventPythia_SLHA)


    // Get next SLHA file path (for use with model CB_SLHA_file_model)
    void getNextSLHAFileName(str& result)
    {
      using namespace Pipes::getNextSLHAFileName;

      static unsigned int counter = 0;
      static bool first = true;

      // Why does this have to be part of the loop?
      if (*Loop::iteration == BASE_INIT)
      {      
        if (first)
        {
          if (!runOptions->hasKey("SLHA_filenames"))
          {
            ColliderBit_error().raise(LOCAL_INFO,"Expected YAML file option 'SLHA_filenames' (a list of SLHA filenames) not found.");
          }
          first = false;
        }

        const static std::vector<str> filenames = runOptions->getValue<std::vector<str> >("SLHA_filenames");

        if (counter >= filenames.size())
        {
          piped_invalid_point.request("No more SLHA files. My work is done.");
          result = "";
        }
        else
        {
          result = filenames.at(counter);
        }
        counter++;          
      }

    }

  }
}
