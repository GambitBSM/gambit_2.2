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
        const str& filename = filenames.at(counter);
        result = std::make_pair(filename, read_SLHA(filename));
      }
      counter++;          

    }


    // Extract SLHA file elements (for use with model CB_SLHA_file_model)
    void getSLHAFileElements(map_str_dbl& result)
    {
      using namespace Pipes::getSLHAFileElements;

      // Split the required SLHAFileNameAndContent pair
      const str& filename = Dep::SLHAFileNameAndContent->first;
      const SLHAstruct& content = Dep::SLHAFileNameAndContent->second;

      // Should missing elements be replaced by a default value?
      const static bool use_missing_element_value = runOptions->hasKey("value_for_missing_elements");
      static double missing_element_value;

      static bool first = true;
      if (first)
      {
        // Check that the required YAML option "SLHA_keys" is present
        if (!runOptions->hasKey("SLHA_keys")) ColliderBit_error().raise(LOCAL_INFO,"Expected YAML file option 'SLHA_keys' (a list of strings, e.g. '- MASS_1000022') not found.");

        // Read default value for missing elements;
        if (use_missing_element_value) missing_element_value = runOptions->getValue<double>("value_for_missing_elements");

        first = false;
      }

      // Read the list of SLHA element keys
      const static std::vector<str> slha_element_keys = runOptions->getValue<std::vector<str> >("SLHA_keys");

      // Loop through the list of SLHA keys and grab the corresponding elements from the SLHA content
      for(str key_str : slha_element_keys)
      {

        // Construct a SLHAea::Key from the key string
        const SLHAea::Key key(key_str);

        // Grab the correct entryand store in the results map
        try
        {
          result[key_str] = SLHAea::to<double>( content.field(key) );
        }
        catch (const std::out_of_range& e)
        {
          std::stringstream errmsg_ss;
          errmsg_ss << "Could not find SLHA element " << key_str << " in file " << filename;

          if (use_missing_element_value)
          {
            logger() << errmsg_ss.str() << EOM;
            result[key_str] = missing_element_value;            
          }
          else
          {
            ColliderBit_error().raise(LOCAL_INFO, errmsg_ss.str());
          }
        }
      }

    }


  }
}
