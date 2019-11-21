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
///  \author Ben Farmer
///          (benjamin.farmer@fysik.su.se)
///    \date 2014 Sep - Dec, 2015 Jan - Mar
///
///  \author Christopher Rogan
///          (christophersrogan@gmail.com)
///  \date 2015 Apr
///
///  \author Tomas Gonzalo
///          (tomas.gonzalo@monash.edu)
///  \date 2016 June, 2019 Oct-Nov
///
///  \author Pat Scott
///          (p.scott@imperial.ac.uk)
///  \date 2015, 2016
///
///  *********************************************

#include "gambit/Elements/gambit_module_headers.hpp"
#include "gambit/Elements/spectrum.hpp"

#include "gambit/SpecBit/SpecBit_rollcall.hpp"
//#include "gambit/SpecBit/SpecBit_helpers.hpp"
#include "gambit/SpecBit/RegisteredSpectra.hpp"

// Switch for debug mode
//#define SPECBIT_DEBUG

namespace Gambit
{

  namespace SpecBit
  {

  /// Check that the spectrum has the canonical LSP for the model being scanned.
    void check_LSP(const Spectrum& spec, const safe_ptr<Options>& runOptions, bool gravitino_is_canonical_LSP)
    {
      double msqd  = spec.get(Par::Pole_Mass, 1000001, 0);
      double msqu  = spec.get(Par::Pole_Mass, 1000002, 0);
      double msl   = spec.get(Par::Pole_Mass, 1000011, 0);
      double msneu = spec.get(Par::Pole_Mass, 1000012, 0);
      double mglui = spec.get(Par::Pole_Mass, 1000021, 0);
      double mchi0 = std::abs(spec.get(Par::Pole_Mass, 1000022, 0));
      double mchip = std::abs(spec.get(Par::Pole_Mass, 1000024, 0));
      double m_canonical_LSP;

      if (gravitino_is_canonical_LSP)
      {
        if (not runOptions->getValueOrDef<bool>(true, "only_gravitino_LSP")) return;
        m_canonical_LSP = spec.get(Par::Pole_Mass, 1000039, 0);
      }
      else
      {
        if (not runOptions->getValueOrDef<bool>(true, "only_neutralino_LSP")) return;
        m_canonical_LSP = mchi0;
      }

      // Check if the canonical LSP in the MSSM version scanned is actually the LSP for this point.
      if (m_canonical_LSP > mchip ||
          m_canonical_LSP > mglui ||
          m_canonical_LSP > msl   ||
          m_canonical_LSP > msneu ||
          m_canonical_LSP > msqu  ||
          m_canonical_LSP > msqd  ||
          m_canonical_LSP > mchi0   )
       {
         str canonical_LSP = (gravitino_is_canonical_LSP ? "Gravitino" : "Neutralino");
         invalid_point().raise(canonical_LSP + " is not LSP.");
       }
    }

    /// Helper to work with pointer
    void check_LSP(const Spectrum* spec, const safe_ptr<Options>& runOptions, bool gravitino_model)
    {
      check_LSP(*spec, runOptions, gravitino_model);
    }

  

    void get_MSSM_spectrum_SPheno (Spectrum& spectrum)
    {
      namespace myPipe = Pipes::get_MSSM_spectrum_SPheno;
      const SMInputs &sminputs = *myPipe::Dep::SMINPUTS;

      // Set up the input structure
      SpectrumInputs inputs(sminputs, SpectrumContents::MSSM(), myPipe::Param, myPipe::runOptions);

      // Get the spectrum from the Backend
      myPipe::BEreq::SPheno_MSSM_Spectrum(spectrum, inputs);

      /// TODO: add LSP check, gravitino mass etc      

      // Add the gravitino mass if it is present
      // TODO: probably this should be an override of sorts
      //if (myPipe::ModelInUse("MSSM63atQ_lightgravitino") or myPipe::ModelInUse("MSSM63atMGUT_lightgravitino"))
      //add_gravitino_mass(slha, myPipe::Param, myPipe::runOptions);

      // Drop SLHA files if requested
      spectrum.drop_SLHAs_if_requested(myPipe::runOptions, "GAMBIT_unimproved_spectrum");

      // Retrieve any mass cuts
      spectrum.check_mass_cuts(*myPipe::runOptions);

    }

    // Runs FlexibleSUSY MSSM spectrum generator with CMSSM (GUT scale) boundary conditions
    // In principle an identical spectrum can be obtained from the function
    // get_MSSMatGUT_spectrum_FS
    // by setting the input parameters to match the CMSSM assumptions
    void get_CMSSM_spectrum_FS (Spectrum& spectrum)
    {
      // Access the pipes for this function to get model and parameter information
      namespace myPipe = Pipes::get_CMSSM_spectrum_FS;

      // Get SLHA2 SMINPUTS values
      const SMInputs& sminputs = *myPipe::Dep::SMINPUTS;

      // Set up the input structure
      SpectrumInputs inputs(sminputs, SpectrumContents::MSSM(), myPipe::Param, myPipe::runOptions);

      // Get the spectrum from the Backend
      myPipe::BEreq::FS_CMSSM_Spectrum(spectrum, inputs);

      /// TODO: add LSP check, gravitino mass etc      

      // Drop SLHA files if requested
      spectrum.drop_SLHAs_if_requested(myPipe::runOptions, "GAMBIT_unimproved_spectrum");

      // Retrieve any mass cuts
      spectrum.check_mass_cuts(*myPipe::runOptions);
    }

    // Runs FlexibleSUSY MSSM spectrum generator with GUT scale input (boundary conditions)
    void get_MSSMatMGUT_spectrum_FS (Spectrum& result)
    {
      // Access the pipes for this function to get model and parameter information
      namespace myPipe = Pipes::get_MSSMatMGUT_spectrum_FS;

      // Get SLHA2 SMINPUTS values
      const SMInputs& sminputs = *myPipe::Dep::SMINPUTS;

      // Set up the input structure
      SpectrumInputs inputs(sminputs, SpectrumContents::MSSM(), myPipe::Param, myPipe::runOptions);

      // TODO: This all goes to the backend source file
      //MSSMatMGUT_input_parameters input;
      //input.mHu2IN = *myPipe::Param.at("mHu2");
      //input.mHd2IN = *myPipe::Param.at("mHd2");
      //input.SignMu = *myPipe::Param.at("SignMu");
      //if(input.SignMu!=-1 and input.SignMu!=1)
      //{
      //   std::ostringstream msg;
      //   msg << "Tried to set SignMu parameter to a value that is not a sign! ("<<input.SignMu<<")! This parameter must be set to either 1 or -1. Please check your inifile and try again.";
      //   SpecBit_error().raise(LOCAL_INFO,msg.str());
      //}
      //fill_MSSM63_input(input,myPipe::Param);

      // Get the spectrum from the Backend
      myPipe::BEreq::FS_MSSMatMGUT_Spectrum(result, inputs);

      /// TODO: add LSP check, gravitino mass etc      

      // Drop SLHA files if requested
      result.drop_SLHAs_if_requested(myPipe::runOptions, "GAMBIT_unimproved_spectrum");

      // Retrieve any mass cuts
      result.check_mass_cuts(*myPipe::runOptions);
    }

    // Runs FlexibleSUSY MSSM spectrum generator with GUT scale input (boundary conditions)
    // but with mA and mu as parameters instead of mHu2 and mHd2
    void get_MSSMatMGUT_mA_spectrum_FS (Spectrum& result)
    {
      // Access the pipes for this function to get model and parameter information
      namespace myPipe = Pipes::get_MSSMatMGUT_mA_spectrum_FS;

      // Get SLHA2 SMINPUTS values
      const SMInputs& sminputs = *myPipe::Dep::SMINPUTS;

      // Set up the input structure
      SpectrumInputs inputs(sminputs, SpectrumContents::MSSM(), myPipe::Param, myPipe::runOptions);

      // TODO: This all goes to the backend source file
      //MSSMatMGUT_mAmu_input_parameters input;
      //input.MuInput  = *myPipe::Param.at("mu");
      // Note this spectrum generator mA2 is the parameter.
      // However this freedom is not used in GAMBIT
      // and not needed as mA is a DRbar mass eigenstate for a scalar
      //double mA = *myPipe::Param.at("mA");
      //input.mA2Input = mA*mA;
      //fill_MSSM63_input(input,myPipe::Param); // Fill the rest

      // Get the spectrum from the Backend
      myPipe::BEreq::FS_MSSMatMGUT_mA_Spectrum(result, inputs);

      /// TODO: add LSP check, gravitino mass etc      

      // Drop SLHA files if requested
      result.drop_SLHAs_if_requested(myPipe::runOptions, "GAMBIT_unimproved_spectrum");

      // Retrieve any mass cuts
      result.check_mass_cuts(*myPipe::runOptions);
    }
 
    // Runs FlexibleSUSY MSSM spectrum generator with EWSB scale input (boundary conditions)
    void get_MSSMatQ_spectrum_FS (Spectrum& result)
    {
      // Access the pipes for this function to get model and parameter information
      namespace myPipe = Pipes::get_MSSMatQ_spectrum_FS;

      // Get SLHA2 SMINPUTS values
      const SMInputs& sminputs = *myPipe::Dep::SMINPUTS;

      // Set up the input structure
      SpectrumInputs inputs(sminputs, SpectrumContents::MSSM(), myPipe::Param, myPipe::runOptions);

      // TODO: This all goes to the backend source file
      //MSSM_input_parameters input;
      //input.Qin    = *myPipe::Param.at("Qin"); // MSSMatQ also requires input scale to be supplied
      //input.mHu2IN = *myPipe::Param.at("mHu2");
      //input.mHd2IN = *myPipe::Param.at("mHd2");
      //input.SignMu = *myPipe::Param.at("SignMu");
      //if(input.SignMu!=-1 and input.SignMu!=1)
      //{
      //   std::ostringstream msg;
      //   msg << "Tried to set SignMu parameter to a value that is not a sign! ("<<input.SignMu<<")! This parameter must be set to either 1 or -1. Please check your inifile and try again.";
      //   SpecBit_error().raise(LOCAL_INFO,msg.str());
      //}
      //fill_MSSM63_input(input,myPipe::Param);

      // Get the spectrum from the Backend
      myPipe::BEreq::FS_MSSMatQ_Spectrum(result, inputs);

      /// TODO: add LSP check, gravitino mass etc      

      // Drop SLHA files if requested
      result.drop_SLHAs_if_requested(myPipe::runOptions, "GAMBIT_unimproved_spectrum");

      // Retrieve any mass cuts
      result.check_mass_cuts(*myPipe::runOptions);
    }

    // Runs FlexibleSUSY MSSM spectrum generator with EWSB scale input (boundary conditions)
    // but with mA and mu as parameters instead of mHu2 and mHd2
    void get_MSSMatQ_mA_spectrum_FS (Spectrum& result)
    {
      // Access the pipes for this function to get model and parameter information
      namespace myPipe = Pipes::get_MSSMatQ_mA_spectrum_FS;

      // Get SLHA2 SMINPUTS values
      const SMInputs& sminputs = *myPipe::Dep::SMINPUTS;

      // Set up the input structure
      SpectrumInputs inputs(sminputs, SpectrumContents::MSSM(), myPipe::Param, myPipe::runOptions);

      // TODO: This all goes to the backend source file
      //MSSM_mAmu_input_parameters input;
      //input.Qin      = *myPipe::Param.at("Qin"); // MSSMatQ also requires input scale to be supplied
      //input.MuInput  = *myPipe::Param.at("mu");
      // Note this spectrum generator mA2 is the parameter.
      // However this freedom is not used in GAMBIT
      // and not needed as mA is a DRbar mass eigenstate for a scalar
      //double mA = *myPipe::Param.at("mA");
      //input.mA2Input = mA*mA;
      //fill_MSSM63_input(input,myPipe::Param); // Fill the rest

      // Get the spectrum from the Backend
      myPipe::BEreq::FS_MSSMatQ_mA_Spectrum(result, inputs);

      /// TODO: add LSP check, gravitino mass etc      

      // Drop SLHA files if requested
      result.drop_SLHAs_if_requested(myPipe::runOptions, "GAMBIT_unimproved_spectrum");

      // Retrieve any mass cuts
      result.check_mass_cuts(*myPipe::runOptions);
    }

    // Runs FlexibleSUSY MSSM spectrum generator with SUSY scale input (boundary conditions)
    // but with mA and mu as parameters instead of mHu2 and mHd2
    void get_MSSMatMSUSY_mA_spectrum_FS (Spectrum& result)
    {
      // Access the pipes for this function to get model and parameter information
      namespace myPipe = Pipes::get_MSSMatMSUSY_mA_spectrum_FS;

      // Get SLHA2 SMINPUTS values
      const SMInputs& sminputs = *myPipe::Dep::SMINPUTS;

      // Set up the input structure
      SpectrumInputs inputs(sminputs, SpectrumContents::MSSM(), myPipe::Param, myPipe::runOptions);

      // TODO: This all goes to the backend source file
      //MSSMatMSUSY_mAmu_input_parameters input;
      //input.MuInput  = *myPipe::Param.at("mu");
      // Note this spectrum generator mA2 is the parameter.
      // However this freedom is not used in GAMBIT
      // and not needed as mA is a DRbar mass eigenstate for a scalar
      //double mA = *myPipe::Param.at("mA");
      //input.mA2Input = mA*mA;    // FS has mA^2 as the parameter
      //fill_MSSM63_input(input,myPipe::Param); // Fill the rest

      // Get the spectrum from the Backend
      myPipe::BEreq::FS_MSSMatMSUSY_mA_Spectrum(result, inputs);

      /// TODO: add LSP check, gravitino mass etc      

      // Drop SLHA files if requested
      result.drop_SLHAs_if_requested(myPipe::runOptions, "GAMBIT_unimproved_spectrum");

      // Retrieve any mass cuts
      result.check_mass_cuts(*myPipe::runOptions);
    }

    // Runs FlexibleSUSY MSSMatMGUTEFTHiggs model spectrum generator
    // and has m3^2 and mu as EWSB outputs, so it is for the
    // MSSMatMGUT_model.
    void get_MSSMatMGUT_spectrum_FlexibleEFTHiggs (Spectrum& result)
    {
      // Access the pipes for this function to get model and parameter information
      namespace myPipe = Pipes::get_MSSMatMGUT_spectrum_FlexibleEFTHiggs;

      // Get SLHA2 SMINPUTS values
      const SMInputs& sminputs = *myPipe::Dep::SMINPUTS;

      // Set up the input structure
      SpectrumInputs inputs(sminputs, SpectrumContents::MSSM(), myPipe::Param, myPipe::runOptions);

      // TODO: This all goes to the backend source file
      //MSSMatMGUTEFTHiggs_input_parameters input;
      //input.mHu2IN = *myPipe::Param.at("mHu2");
      //input.mHd2IN = *myPipe::Param.at("mHd2");
      //input.SignMu = *myPipe::Param.at("SignMu");
      //if(input.SignMu!=-1 and input.SignMu!=1)
      //{
      //   std::ostringstream msg;
      //   msg << "Tried to set SignMu parameter to a value that is not a sign! ("<<input.SignMu<<")! This parameter must be set to either 1 or -1. Please check your inifile and try again.";
      //   SpecBit_error().raise(LOCAL_INFO,msg.str());
      //}

      //fill_MSSM63_input(input,myPipe::Param); // Fill the rest

      // Get the spectrum from the Backend
      myPipe::BEreq::FS_MSSMatMGUTEFTHiggs_Spectrum(result, inputs);

      /// TODO: add LSP check, gravitino mass etc      

      // Drop SLHA files if requested
      result.drop_SLHAs_if_requested(myPipe::runOptions, "GAMBIT_unimproved_spectrum");

      // Retrieve any mass cuts
      result.check_mass_cuts(*myPipe::runOptions);
    }

    // Runs FlexibleSUSY MSSMatMGUTEFTHiggs model spectrum generator
    // and has m3^2 and mu as EWSB outputs, so it is for the
    // MSSMatMGUT_model.
    void get_MSSMatMGUT_mA_spectrum_FlexibleEFTHiggs (Spectrum& result)
    {
      // Access the pipes for this function to get model and parameter information
      namespace myPipe = Pipes::get_MSSMatMGUT_mA_spectrum_FlexibleEFTHiggs;

      // Get SLHA2 SMINPUTS values
      const SMInputs& sminputs = *myPipe::Dep::SMINPUTS;

      // Set up the input structure
      SpectrumInputs inputs(sminputs, SpectrumContents::MSSM(), myPipe::Param, myPipe::runOptions);

      // TODO: This all goes to the backend source file
      //MSSMatMGUTEFTHiggs_mAmu_input_parameters input;
      //input.MuInput  = *myPipe::Param.at("mu");
      // Note this spectrum generator mA2 is the parameter.
      // However this freedom is not used in GAMBIT
      // and not needed as mA is a DRbar mass eigenstate for a scalar
      //double mA = *myPipe::Param.at("mA");
      //input.mA2Input = mA*mA;

      //fill_MSSM63_input(input,myPipe::Param); // Fill the rest

      // Get the spectrum from the Backend
      myPipe::BEreq::FS_MSSMatMGUTEFTHiggs_mA_Spectrum(result, inputs);

      /// TODO: add LSP check, gravitino mass etc      

      // Drop SLHA files if requested
      result.drop_SLHAs_if_requested(myPipe::runOptions, "GAMBIT_unimproved_spectrum");

      // Retrieve any mass cuts
      result.check_mass_cuts(*myPipe::runOptions);
    }

    // Runs FlexibleSUSY MSSMEFTHiggs model spectrum generator
    // and has m3^2 and mu as EWSB outputs, so it is for the
    // MSSMatQ_model.
    void get_MSSMatQ_spectrum_FlexibleEFTHiggs (Spectrum& result)
    {
      // Access the pipes for this function to get model and parameter information
      namespace myPipe = Pipes::get_MSSMatQ_spectrum_FlexibleEFTHiggs;

      // Get SLHA2 SMINPUTS values
      const SMInputs& sminputs = *myPipe::Dep::SMINPUTS;

      // Set up the input structure
      SpectrumInputs inputs(sminputs, SpectrumContents::MSSM(), myPipe::Param, myPipe::runOptions);

      // TODO: This all goes to the backend source file
      //MSSMEFTHiggs_input_parameters input;
      // MSSMatQ also requires input scale to be supplied with name MSUSY
      //input.MSUSY  = *myPipe::Param.at("Qin");
      //input.mHu2IN = *myPipe::Param.at("mHu2");
      //input.mHd2IN = *myPipe::Param.at("mHd2");
      //input.SignMu = *myPipe::Param.at("SignMu");
      //fill_MSSM63_input_altnames(input,myPipe::Param); // Fill the rest

      // Get the spectrum from the Backend
      myPipe::BEreq::FS_MSSMatQEFTHiggs_Spectrum(result, inputs);

      /// TODO: add LSP check, gravitino mass etc      

      // Drop SLHA files if requested
      result.drop_SLHAs_if_requested(myPipe::runOptions, "GAMBIT_unimproved_spectrum");

      // Retrieve any mass cuts
      result.check_mass_cuts(*myPipe::runOptions);
    }

    // Runs FlexibleSUSY MSSMEFTHiggs_mAmu spectrum generator with
    // boundary conditions at a user specified scale, ie accepts MSSM
    // parameters at Q, and has DRbar mA and mu as an input and mHu2
    // and mHd2 as EWSB outputs, so it is for the MSSMatMSUSY_mA model.
    void get_MSSMatQ_mA_spectrum_FlexibleEFTHiggs (Spectrum& result)
    {
      // Access the pipes for this function to get model and parameter information
      namespace myPipe = Pipes::get_MSSMatQ_mA_spectrum_FlexibleEFTHiggs;

      // Get SLHA2 SMINPUTS values
      const SMInputs& sminputs = *myPipe::Dep::SMINPUTS;

      // Set up the input structure
      SpectrumInputs inputs(sminputs, SpectrumContents::MSSM(), myPipe::Param, myPipe::runOptions);

      // TODO: This all goes to the backend source file
      //MSSMEFTHiggs_mAmu_input_parameters input;
      //input.MuInput  = *myPipe::Param.at("mu");
      // This FS spectrum generator has mA as the parameter
      //input.mAInput = *myPipe::Param.at("mA");
      // Note: Qin has been named MSUSY inside the spectrum generator
      // but it is a user-input scale in this case.
      //input.MSUSY = *myPipe::Param.at("Qin");
      // Fill the rest.
      // Note: This particular spectrum generator has been created with
      // different names for parameter inputs.  We should standardise this
      //fill_MSSM63_input_altnames(input,myPipe::Param);

      // Get the spectrum from the Backend
      myPipe::BEreq::FS_MSSMatQEFTHiggs_mA_Spectrum(result, inputs);

      /// TODO: add LSP check, gravitino mass etc      

      // Drop SLHA files if requested
      result.drop_SLHAs_if_requested(myPipe::runOptions, "GAMBIT_unimproved_spectrum");

      // Retrieve any mass cuts
      result.check_mass_cuts(*myPipe::runOptions);
    }
 
    // Runs FlexibleSUSY MSSMEFTHiggs model spectrum generator with SUSY
    // scale boundary conditions, ie accepts MSSM parameters at MSUSY,
    // and has DRbar mA and mu as an input and mHu2 and mHd2 as EWSB
    // outputs, so it is for the MSSMatMSUSY_mA model.
    void get_MSSMatMSUSY_mA_spectrum_FlexibleEFTHiggs (Spectrum& result)
    {
      // Access the pipes for this function to get model and parameter information
      namespace myPipe = Pipes::get_MSSMatMSUSY_mA_spectrum_FlexibleEFTHiggs;

      // Get SLHA2 SMINPUTS values
      const SMInputs& sminputs = *myPipe::Dep::SMINPUTS;

      // Set up the input structure
      SpectrumInputs inputs(sminputs, SpectrumContents::MSSM(), myPipe::Param, myPipe::runOptions);

      // TODO: This all goes to the backend source file
      //MSSMatMSUSYEFTHiggs_mAmu_input_parameters input;
      //input.MuInput  = *myPipe::Param.at("mu");
      // This FS spectrum generator has mA as the parameter
      //input.mAInput = *myPipe::Param.at("mA");
      //fill_MSSM63_input_altnames(input,myPipe::Param); // Fill the rest

      // Get the spectrum from the Backend
      myPipe::BEreq::FS_MSSMatMSUSYEFTHiggs_mA_Spectrum(result, inputs);

      /// TODO: add LSP check, gravitino mass etc      

      // Drop SLHA files if requested
      result.drop_SLHAs_if_requested(myPipe::runOptions, "GAMBIT_unimproved_spectrum");

      // Retrieve any mass cuts
      result.check_mass_cuts(*myPipe::runOptions);
    }

  
  } // end namespace SpecBit
} // end namespace Gambit

