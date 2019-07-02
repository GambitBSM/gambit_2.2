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


    // Get next SLHA file path and content (for use with model CB_SLHA_file_model)
    void getNextSLHAFileNameAndContent(pair_str_SLHAstruct& result)
    {
      using namespace Pipes::getNextSLHAFileNameAndContent;

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
          result = std::make_pair("", SLHAstruct());
        }
        else
        {

          // logger() << "Reading SLHA file: " << *Dep::SLHAFileName << EOM;               \
          // std::ifstream ifs(Dep::SLHAFileName->c_str());                                \
          // if(!ifs.good()){ ColliderBit_error().raise(LOCAL_INFO,"ERROR: SLHA file not found."); } \
          // ifs >> slha;                                                                  \

          const str& filename = filenames.at(counter);
          result = std::make_pair(filename, read_SLHA(filename));
        }
        counter++;          
      }

    }


   /*
  
    // Extract SLHA file elements (for use with model CB_SLHA_file_model)
    void getSLHAFileElements(map_str_dbl& result)
    {
      using namespace getSLHAFileElements;

      static int counter = 0
      counter++;

      const str& filename = Dep::SLHAFileNameAndContent->first;
      const SLHAstruct& content = Dep::SLHAFileNameAndContent->second;

      cout << "DEBUG: getSLHAFileElements: counter: " << counter << " filename: " << filename << endl;

      // Check that the required YAML option is present
      static bool first = true;
      if (first)
      {
        if (!runOptions->hasKey("SLHA_elements"))
        {
          ColliderBit_error().raise(LOCAL_INFO,"Expected YAML file option 'SLHA_elements' (a list of strings, e.g. '- MASS_1000022') not found.");
        }
        first = false;
      }

      // Read the list of SLHA element keys
      const static std::vector<str> slha_element_keys = runOptions->getValue<std::vector<str> >("SLHA_elements");

      for(str key : slha_element_keys)
      {
        key = "NMIX_1_1"

        std::replace(key.begin(), key.end(), ':', ' ');        
        vector<str> temp_str_vector;
        vector<str> temp_str_vector;
        std::stringstream ss(key);
        while (ss >> )

  string str = "102:330:3133:76531:451:000:12:44412";
std::replace(str.begin(), str.end(), ':', ' ');  // replace ':' by ' '

vector<int> array;
stringstream ss(str);
int temp;
while (ss >> temp)
    array.push_back(temp); 

      }

        // For each key:
        //   - construct two variables:
        //       - block_name : "MASS", "NMIX", "DECAY"
        //       - indices : [1000024], [1,1], [1000024, 1000022, 24] 
        //
        //   - check is_decay and length of index vector
        //
        //   - get the element value from the SLHAstruct (how to get the DECAY...?) 
        //
        //   - store in result map: {
        //                           "MASS_1000024": 134.23, 
        //                           "NMIX_1_1": 0.32,
        //                           "DECAY_1000024_1000022_24": 1.0
        //                          }
        //
        // Use static variables to avoid having to decompose the key strings for every point
        //


    }

  */

  }
}
