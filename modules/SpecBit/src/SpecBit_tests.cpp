//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  Functions of module SpecBit
///
///  These functions link ModelParameters to 
///  Spectrum objects in various ways.
///
///  *********************************************
///
///  Authors (add name and date if you modify):
///
///  \author Ben Farmer
///          (ben.farmer@gmail.com)
///    \date 2014 Sep
///  
///  *********************************************

#include "gambit_module_headers.hpp"
#include "SpecBit_rollcall.hpp"

// Flexible SUSY stuff (should not be needed by the rest of gambit)
#include "CMSSM_two_scale_model.hpp"
#include "CMSSM_two_scale_model_slha.hpp"
//#include "CMSSM_physical.hpp"

#include "ew_input.hpp"
//#include "logger.hpp"
//#include "wrappers.hpp"
#include "MSSMSpec.hpp"
#include "numerics.hpp"

// Used in SpecBit_tests.hpp to switch the test output to logger() rather than std::cout when run through SpecBit.
#define IN_SPECBIT
// Spectrum object test functions
#include "SpecBit_tests.hpp"
#include "SpecBit_examples.hpp"
namespace Gambit
{

  namespace SpecBit
  {
    using namespace LogTags;
    using namespace flexiblesusy;

    /// Quick macro to simplify the check of Pipe::Models
    #define QUERYMODELS(MODEL) std::find(Pipe::Models->begin(), Pipe::Models->end(), MODEL) != Pipe::Models->end()
    
    /// Create a spectrum object for testing purposes
    void make_test_spectrum(Spectrum* &result)
    {
      typedef CMSSM_interface<ALGORITHM1> MI;
      static MI::Model FS_model; //start with empty flexiblesusy object
      // Or could use flexiblesusy classes directly; these two are equivalent in this case:
      //static CMSSM_slha<Two_scale> FS_model; //start with empty flexiblesusy object

      // Create model interface class (leaving input stuff with default values)
      MI model_interface(FS_model);

      // Create Spectrum object to wrap flexiblesusy object
      static MSSMSpec<MI> mssm(FS_model);

      // I think these objects should only get created once since they are static...      
      // ...and they should be destructed automatically when the program ends.

      setup(mssm.model); //fill with some parameters
      mssm.model.calculate_DRbar_parameters(); //calculated DRbar masses 
      mssm.model.calculate_pole_masses();//now calculate pole masses

      // Check contents
      logger() << "This is specbit_tests. Checking Spectrum object contents..." << std::endl;
      if(TestMssmParGets(mssm, mssm.model)==false){
          logger() << "TestMssmParGets fail." << std::endl;
          return;
       }
       if(TestMssmPoleGets(mssm, mssm.model)==false){
          logger() << "TestMssmPoleGets fail." << std::endl;
          return;
       }
       //So now we have a mssm1 model object filled, as it will be
       //stored in Gambit after the spectrum generator has run
       // mssm.mass2_par_mapping(); //call mapping - this needs to be changed.
    
       mssm_manipulate(mssm);  //function can manipulate knowing the model
 
      // Store result for gambit to use
      result = &mssm;
    }

    /// Function to test out SpecBit features
    void specbit_test_func1 (double &result)
    {
      // Access the pipes for this function to get model and parameter information
      using namespace Pipes::specbit_test_func1;

      std::cout << "Running specbit_test_func1" << std::endl;

      Spectrum* spec = *Dep::MSSM_spectrum; //Test retrieve pointer to Spectrum object 

      //const Spectrum& spec(*(Dep::particle_spectrum->get())); // Get Spectrum object ptr out of dependency pipe and make a nice reference out of it.
      spec_manipulate(spec); //function can manipulate without knowing model.

      logger() << EOM;  
    }

    /// Function to test out SpecBit features
    void specbit_test_func2 (double &result)
    {
      std::cout << "Running specbit_test_func2" << std::endl;

      // TESTING
      // Direct access to flexiblesusy function, for testing
      CMSSM_slha<Two_scale> FS_model; //start with empty flexiblesusy object

      // Create model interface class (leaving input stuff with default values)
      CMSSM_interface<Two_scale> model_interface(FS_model);

      // Create Spectrum object to wrap flexiblesusy object
      MSSMSpec<CMSSM_interface<Two_scale>> mssm(model_interface);

      // Test run functions
      std::cout << "Spectrum via MSSMSpec" << std::endl;
      std::cout << "mssm.runningpars.GetScale() =" 
          << mssm.runningpars.GetScale() << std::endl;
      std::cout << "mHd2 = "  
          << mssm.runningpars.get_mass2_parameter("mHd2") << std::endl;
      std::cout << "mHu2 = "  
          << mssm.runningpars.get_mass2_parameter("mHu2") << std::endl;

      // Do it again using a Spectrum base pointer
      Spectrum* spec = &mssm;
      std::cout << "Spectrum via Spectrum*" << std::endl;
      std::cout << "spec->runningpars.GetScale() =" 
          << spec->runningpars.GetScale() << std::endl;
      std::cout << "mHd2 = "  
          << spec->runningpars.get_mass2_parameter("mHd2") << std::endl;
      std::cout << "mHu2 = "  
          << spec->runningpars.get_mass2_parameter("mHu2") << std::endl;

      // Fill the model and do it again
      std::cout << "Spectrum via Spectrum* (filled)" << std::endl;
      setup(mssm.model);
      std::cout << "spec->runningpars.GetScale() =" 
          << spec->runningpars.GetScale() << std::endl;
      std::cout << "mHd2 = "  
          << spec->runningpars.get_mass2_parameter("mHd2") << std::endl;
      std::cout << "mHu2 = "  
          << spec->runningpars.get_mass2_parameter("mHu2") << std::endl;

    }

    /// Function to test out SpecBit features
    void specbit_test_func3 (double &result)
    {
      // Requests a Spectrum object of capability SM_spectrum; test what we can retrieve from this
      using namespace Pipes::specbit_test_func3;
      Spectrum* spec = *Dep::SM_spectrum; //Test retrieve pointer to Spectrum object 

      SM_checks(spec); // Run some tests on standard model parameters 
      logger() << EOM;
    }


  } // end namespace SpecBit
} // end namespace Gambit

