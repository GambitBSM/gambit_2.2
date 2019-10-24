//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  Functions of module SpecBit, model MSSM
///
///  *********************************************
///
///  Authors (add name and date if you modify):
///
///  \author Tomas Gonzalo
///          (tomas.gonzalo@monash.edu)
///  \date 2019
///
///  *********************************************

#include "gambit/Elements/gambit_module_headers.hpp"
//#include "gambit/Elements/spectrum_factories.hpp"
//#include "gambit/Elements/smlike_higgs.hpp"
//#include "gambit/Models/SimpleSpectra/MSSMSimpleSpec.hpp"
//#include "gambit/Utils/stream_overloads.hpp" // Just for more convenient output to logger
//#include "gambit/Utils/util_macros.hpp"
//#include "gambit/Utils/util_types.hpp" // struct SpectrumInputs
#include "gambit/SpecBit/SpecBit_rollcall.hpp"
//#include "gambit/SpecBit/SpecBit_helpers.hpp"
//#include "gambit/SpecBit/QedQcdWrapper.hpp"
//#include "gambit/SpecBit/MSSMSpec.hpp"
//#include "gambit/SpecBit/model_files_and_boxes.hpp" // #includes lots of flexiblesusy headers and defines interface classes

// Switch for debug mode
//#define SPECBIT_DEBUG

namespace Gambit
{

  namespace SpecBit
  {

    void get_MSSM_spectrum_SPheno (Spectrum& spectrum)
    {
      namespace myPipe = Pipes::get_MSSM_spectrum_SPheno;
      const SMInputs &sminputs = *myPipe::Dep::SMINPUTS;

      // Set up the input structure
      SpectrumInputs inputs;
      inputs.sminputs = sminputs;
      inputs.param = myPipe::Param;
      inputs.options = myPipe::runOptions;

      // Retrieve any mass cuts
      /// TODO: These should be a struct in spectrum class
      //static const Spectrum::mc_info mass_cut = myPipe::runOptions->getValueOrDef<Spectrum::mc_info>(Spectrum::mc_info(), "mass_cut");
      //static const Spectrum::mr_info mass_ratio_cut = myPipe::runOptions->getValueOrDef<Spectrum::mr_info>(Spectrum::mr_info(), "mass_ratio_cut");

      // Get the spectrum from the Backend
      myPipe::BEreq::SPheno_MSSM_Spectrum(spectrum, inputs);

      // Add the gravitino mass if it is present
      // TODO: probably this should be an override of sorts
      //if (myPipe::ModelInUse("MSSM63atQ_lightgravitino") or myPipe::ModelInUse("MSSM63atMGUT_lightgravitino"))
      //add_gravitino_mass(slha, myPipe::Param, myPipe::runOptions);

    }

    // Runs FlexibleSUSY MSSM spectrum generator with CMSSM (GUT scale) boundary conditions
    void get_CMSSM_spectrum_FS (Spectrum& spectrum)
    {

      // Access the pipes for this function to get model and parameter information
      namespace myPipe = Pipes::get_CMSSM_spectrum_FS;

      // Get SLHA2 SMINPUTS values
      const SMInputs& sminputs = *myPipe::Dep::SMINPUTS;

      // Set up the input structure
      SpectrumInputs inputs;
      inputs.sminputs = sminputs;
      inputs.param = myPipe::Param;
      inputs.options = myPipe::runOptions;

      // Get the spectrum from the Backend
      myPipe::BEreq::FS_CMSSM_Spectrum(spectrum, inputs);

      // Drop SLHA files if requested
      spectrum.drop_SLHAs_if_requested(myPipe::runOptions, "GAMBIT_unimproved_spectrum");

    }

  } // end namespace SpecBit
} // end namespace Gambit

