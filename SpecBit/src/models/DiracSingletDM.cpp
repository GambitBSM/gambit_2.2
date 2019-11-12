//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  Functions of module SpecBit
///
///  SpecBit module functions related to the
///  DiracSingletDM model
///
///  *********************************************
///
///  Authors (add name and date if you modify):
///
///  \author Ankit Beniwal
///          (ankit.beniwal@adelaide.edu.au)
///  \date 2016 Oct, Nov
///  \date 2017 Jun, Sep
///  \date 2018 Feb
///
///  \author Sanjay Bloor
///          (sanjay.bloor12@imperial.ac.uk)
///  \date 2018 Aug
///
///  \author Tomas Gonzalo
///          (tomas.gonzalo@monash.edu)
///  \date 2019 Nov
///
///  *********************************************

#include <string>
#include <sstream>

#include "gambit/Elements/gambit_module_headers.hpp"

#include "gambit/Elements/spectrum.hpp"

#include "gambit/SpecBit/SpecBit_rollcall.hpp"
#include "gambit/SpecBit/RegisteredSpectra.hpp"
#include "gambit/Models/partmap.hpp"

// Switch for debug mode
//#define SPECBIT_DEBUG

namespace Gambit
{

  namespace SpecBit
  {
    using namespace LogTags;

    /// Get a (simple) Spectrum object wrapper for the DiracSingletDM_Z2 model
    void get_DiracSingletDM_Z2_spectrum(Spectrum& result)
    {
      namespace myPipe = Pipes::get_DiracSingletDM_Z2_spectrum;
      const SMInputs& sminputs = *myPipe::Dep::SMINPUTS;

      // Initialise an SLHAea object to carry the Dirac plus Higgs sector information
      SLHAea::Coll slha;
      SLHAea_add_block(slha, "MASS");
      SLHAea_add_block(slha, "VEVS");
      SLHAea_add_block(slha, "COUPLINGS");
      SLHAea_add_block(slha, "SINTHETAW");
      SLHAea_add_block(slha, "GAUGE");
      SLHAea_add_block(slha, "YD");
      SLHAea_add_block(slha, "YU");
      SLHAea_add_block(slha, "YE");

      // quantities needed to fill container spectrum, intermediate calculations
      double alpha_em = 1.0 / sminputs.alphainv;
      double C = alpha_em * pi / (sminputs.GF * pow(2,0.5));
      double sinW2 = 0.5 - pow( 0.25 - C/pow(sminputs.mZ,2) , 0.5);
      double cosW2 = 0.5 + pow( 0.25 - C/pow(sminputs.mZ,2) , 0.5);
      double e = pow( 4*pi*( alpha_em ),0.5) ;

      // Higgs sector
      double mh   = *myPipe::Param.at("mH");
      std::pair<int,int> pdg = Models::ParticleDB().pdg_pair("h0_1");
      slha["MASS"][""] << pdg.first << mh << "# mH";

      double vev        = 1. / sqrt(sqrt(2.)*sminputs.GF);
      slha["VEVS"][""] << 1 << vev << " # vev";
      // slha["COUPLINGS"][""] << 2 << GF*pow(mh,2)/pow(2,0.5) << " # lambda_h";

      // DiracSingletDM_Z2 sector
      pdg = Models::ParticleDB().pdg_pair("F");
      double mF = *myPipe::Param.at("mF");
      slha["MASS"][""] << pdg.first << mF << "# mF";
      double lF = *myPipe::Param.at("lF");
      slha["COUPLINGS"][""] << 1 << lF << "# lF";
      double xi = *myPipe::Param.at("xi");
      slha["COUPLINGS"][""] << 3 << xi << "# xi";

      // Invalidate point if the EFT validity constraint is not satisfied
      // See https://arxiv.org/abs/1512.06458v4 for more details
      if (myPipe::runOptions->getValueOrDef<bool>(false,"impose_EFT_validity"))
      {
        // Different EFT validity constraints for different model parametrisations.
        if (myPipe::ModelInUse("DiracSingletDM_Z2_sps"))
        {
          // Invadlidate point if the EFT validity constraint is not satisfied,
          // for each coupling independently.
          double gs = lF * std::cos(xi);
          double gp = lF * std::sin(xi);

          if (myPipe::runOptions->getValueOrDef<bool>(false,"impose_EFT_validity"))
          {
            if (gs >= (4*pi)/(2*mF))
            {
              std::ostringstream msg;
              msg << "Parameter point [mF, lF_s] = [" << mF << " GeV, "
                  << gs << "/GeV] does not satisfy the EFT validity constraint.";
              invalid_point().raise(msg.str());
            }
            if (gp >= (4*pi)/(2*mF))
            {
              std::ostringstream msg;
              msg << "Parameter point [mF, lF_ps] = [" << mF << " GeV, "
                  << gp << "/GeV] does not satisfy the EFT validity constraint.";
              invalid_point().raise(msg.str());
            }
          }
        }
        else
        {
          // Parametrisation with lambda/Lambda, xi
          if (lF >= (4*pi)/(2*mF))
          {
            std::ostringstream msg;
            msg << "Parameter point [mF, lF] = [" << mF << " GeV, "
                << lF << "/GeV] does not satisfy the EFT validity constraint.";
            invalid_point().raise(msg.str());
          }
        }
      }

      // Standard model
      slha["SINTHETAW"][""] << 1 << sinW2 << "# sinW2";

      // gauge couplings
      slha["GAUGE"][""] << 1 << e / sqrt(sinW2) << "# g1";
      slha["GAUGE"][""] << 2 << e / sqrt(cosW2) << "# g2";
      slha["GAUGE"][""] << 3 << pow( 4*pi*( sminputs.alphaS ),0.5) << "# g3";

      // Yukawas
      double sqrt2v = pow(2.0,0.5)/vev;
      slha["YU"][""] << 1 << 1 << sqrt2v * sminputs.mU << "# Yu(1,1)";
      slha["YU"][""] << 2 << 2 << sqrt2v * sminputs.mCmC << "# Yu(2,2)";
      slha["YU"][""] << 3 << 3 << sqrt2v * sminputs.mT << "# Yu(3,3)";
      slha["YE"][""] << 1 << 1 << sqrt2v * sminputs.mE << "# Ye(1,1)";
      slha["YE"][""] << 2 << 2 << sqrt2v * sminputs.mMu << "# Ye(2,2)";
      slha["YE"][""] << 3 << 3 << sqrt2v * sminputs.mTau << "# Ye(3,3)";
      slha["YD"][""] << 1 << 1 << sqrt2v * sminputs.mD << "# Yd(1,1)";
      slha["YD"][""] << 2 << 2 << sqrt2v * sminputs.mS << "# Yd(2,2)";
      slha["YD"][""] << 3 << 3 << sqrt2v * sminputs.mBmB << "# Yd(3,3)";

      // SpectrumContents struct
      SpectrumContents::DiracSingletDM_Z2 dirac;

      // Create spectrum object
      // Take mZ as the spectrum scale
      result = Spectrum(slha, dirac, sminputs.mZ, false);

      // Retrieve any mass cuts
      result.check_mass_cuts(*myPipe::runOptions);

    }

    // print spectrum out, stripped down copy from MSSM version with variable names changed
    void fill_map_from_DiracSingletDM_Z2spectrum(std::map<std::string,double>& specmap, const Spectrum& diracdmspec)
    {
      /// Add everything... use spectrum contents routines to automate task
      static const SpectrumContents::DiracSingletDM_Z2 contents;
      static const std::vector<SpectrumContents::Parameter> required_parameters = contents.all_parameters();

      for(std::vector<SpectrumContents::Parameter>::const_iterator it = required_parameters.begin();
           it != required_parameters.end(); ++it)
      {
         const Par::Tags        tag   = it->tag();
         const std::string      name  = it->name();
         const std::vector<int> shape = it->shape();

         /// Verification routine should have taken care of invalid shapes etc, so won't check for that here.

         // Check scalar case
         if(shape.size()==1 and shape[0]==1)
         {
           std::ostringstream label;
           label << name <<" "<< Par::toString.at(tag);
           specmap[label.str()] = diracdmspec.get(tag,name);
         }
         // Check vector case
         else if(shape.size()==1 and shape[0]>1)
         {
           for(int i = 1; i<=shape[0]; ++i) {
             std::ostringstream label;
             label << name <<"_"<<i<<" "<< Par::toString.at(tag);
             specmap[label.str()] = diracdmspec.get(tag,name,i);
           }
         }
         // Check matrix case
         else if(shape.size()==2)
         {
           for(int i = 1; i<=shape[0]; ++i) {
             for(int j = 1; j<=shape[0]; ++j) {
               std::ostringstream label;
               label << name <<"_("<<i<<","<<j<<") "<<Par::toString.at(tag);
               specmap[label.str()] = diracdmspec.get(tag,name,i,j);
             }
           }
         }
         // Deal with all other cases
         else
         {
           // ERROR
           std::ostringstream errmsg;
           errmsg << "Error, invalid parameter received while converting DiracSingletDM_Z2spectrum to map of strings! This should no be possible if the spectrum content verification routines were working correctly; they must be buggy, please report this.";
           errmsg << "Problematic parameter was: "<< tag <<", " << name << ", shape="<< shape;
           utils_error().forced_throw(LOCAL_INFO,errmsg.str());
         }
      }

    }

    void get_DiracSingletDM_Z2_spectrum_as_map (std::map<std::string,double>& specmap)
    {
      namespace myPipe = Pipes::get_DiracSingletDM_Z2_spectrum_as_map;
      const Spectrum& diracdmspec(*myPipe::Dep::DiracSingletDM_Z2_spectrum);
      fill_map_from_DiracSingletDM_Z2spectrum(specmap, diracdmspec);
    }

  } // end namespace SpecBit
} // end namespace Gambit
