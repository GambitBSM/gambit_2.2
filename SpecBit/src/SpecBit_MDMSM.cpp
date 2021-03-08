//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  Implementation of SpecBit routines for 
///  MDMSM.
///
///  Authors (add name and date if you modify):    
///       *** Automatically created by GUM ***     
///                                                
///  \author The GAMBIT Collaboration             
///  \date 04:53PM on March 08, 2021
///                                                
///  ********************************************* 

#include <string>
#include <sstream>

#include "gambit/Elements/gambit_module_headers.hpp"
#include "gambit/Elements/spectrum.hpp"
#include "gambit/Elements/spectrum_factories.hpp"
#include "gambit/Utils/stream_overloads.hpp"
#include "gambit/Utils/util_macros.hpp"

#include "gambit/SpecBit/SpecBit_rollcall.hpp"
#include "gambit/SpecBit/SpecBit_helpers.hpp"
#include "gambit/Models/SpectrumContents/RegisteredSpectra.hpp"
#include "gambit/SpecBit/QedQcdWrapper.hpp"
#include "gambit/Models/SimpleSpectra/MDMSMSimpleSpec.hpp"
#include "gambit/Models/SimpleSpectra/SMHiggsSimpleSpec.hpp"


namespace Gambit
{
  
  namespace SpecBit
  {
    using namespace LogTags;

    /// A function to calculate the effective ggY coupling strength
    std::tuple<double,double> getggY(double mt, double mY, double as, double vev, double cY) {

      double tau = 4 * mt*mt/(mY*mY);
      std::complex<double> ftau;
      std::complex<double> i = sqrt(-1.);

      if (tau < 1.) {
        ftau = -0.25 * pow(log((1+sqrt(1-tau))/(1 - sqrt(1-tau))) - i * M_PI,2);
      } else {
        ftau = pow(asin(1/(sqrt(tau))),2);
      }

      std::complex<double> one(1.0,0.0);
      double gggY1 = (as * cY / (8 * M_PI * vev) * tau * (one + (1.-tau) * ftau)).real(); 
      double gggY2 = (as * cY / (8 * M_PI * vev) * tau * (one + (1.-tau) * ftau)).imag(); 

      return(std::make_tuple(gggY1,gggY2));

    }
    
    /// Get a (simple) Spectrum object wrapper for MDMSM_spectrum model.
    void get_MDMSM_spectrum(Spectrum& result)
    {
      namespace myPipe = Pipes::get_MDMSM_spectrum;
      const SMInputs& sminputs = *myPipe::Dep::SMINPUTS;
      
      // Initialise SLHAea object 
      SLHAstruct slha;
      
      // Block DMINT
      SLHAea_add_block(slha, "DMINT");
      SLHAea_add(slha, "DMINT", 1, *myPipe::Param["gchi"]);
      SLHAea_add(slha, "DMINT", 2, *myPipe::Param["cY"]);
//      SLHAea_add(slha, "DMINT", 3, *myPipe::Param["gggY1"]);
//      SLHAea_add(slha, "DMINT", 4, *myPipe::Param["gggY2"]);
      double vev = 1. / sqrt(sqrt(2.)*sminputs.GF);
      double sqrt2v = pow(2.0,0.5)/vev;
      
      SLHAea_add_block(slha, "VEVS");
      SLHAea_add(slha, "VEVS", 1, vev);
      
      SLHAea_add_block(slha, "HMIX");
      SLHAea_add(slha, "HMIX", 3, vev);
      
      
      // Block MASS
      SLHAea_add_block(slha, "MASS");
      SLHAea_add(slha, "MASS", 25, *myPipe::Param["mH"]);
      SLHAea_add(slha, "MASS", 52, *myPipe::Param["mchi"]);
      SLHAea_add(slha, "MASS", 99902, *myPipe::Param["mY"]);
      
      // quantities needed to fill container spectrum
      double alpha_em = 1.0 / sminputs.alphainv;
      double C = alpha_em * pi / (sminputs.GF * pow(2,0.5));
      double sinW2 = 0.5 - pow( 0.25 - C/pow(sminputs.mZ,2) , 0.5);
      double cosW2 = 0.5 + pow( 0.25 - C/pow(sminputs.mZ,2) , 0.5);
      double e = pow( 4*pi*( alpha_em ),0.5);
      
      SLHAea_add_block(slha, "GAUGE");
      SLHAea_add(slha, "GAUGE", 1, sqrt(5/3) * e / sqrt(cosW2) );
      SLHAea_add(slha, "GAUGE", 2, e / sqrt(sinW2));
      SLHAea_add(slha, "GAUGE", 3, pow( 4*pi*sminputs.alphaS,0.5) );
      
      SLHAea_add_block(slha, "SINTHETAW");
      SLHAea_add(slha, "SINTHETAW", 1, sinW2);
      
      SLHAea_add_block(slha, "YU");
      SLHAea_add(slha, "YU", 1, 1, sqrt2v*sminputs.mU, "u");
      SLHAea_add(slha, "YU", 1, 2, 0., "");
      SLHAea_add(slha, "YU", 1, 3, 0., "");
      SLHAea_add(slha, "YU", 2, 1, 0., "");
      SLHAea_add(slha, "YU", 2, 2, sqrt2v*sminputs.mCmC, "c");
      SLHAea_add(slha, "YU", 2, 3, 0., "");
      SLHAea_add(slha, "YU", 3, 1, 0., "");
      SLHAea_add(slha, "YU", 3, 2, 0., "");
      SLHAea_add(slha, "YU", 3, 3, sqrt2v*sminputs.mT, "t");
      
      SLHAea_add_block(slha, "YE");
      SLHAea_add(slha, "YE", 1, 1, sqrt2v*sminputs.mE, "e");
      SLHAea_add(slha, "YE", 1, 2, 0., "");
      SLHAea_add(slha, "YE", 1, 3, 0., "");
      SLHAea_add(slha, "YE", 2, 1, 0., "");
      SLHAea_add(slha, "YE", 2, 2, sqrt2v*sminputs.mMu, "mu");
      SLHAea_add(slha, "YE", 2, 3, 0., "");
      SLHAea_add(slha, "YE", 3, 1, 0., "");
      SLHAea_add(slha, "YE", 3, 2, 0., "");
      SLHAea_add(slha, "YE", 3, 3, sqrt2v*sminputs.mTau, "tau");
      
      SLHAea_add_block(slha, "YD");
      SLHAea_add(slha, "YD", 1, 1, sqrt2v*sminputs.mD, "d");
      SLHAea_add(slha, "YD", 1, 2, 0., "");
      SLHAea_add(slha, "YD", 1, 3, 0., "");
      SLHAea_add(slha, "YD", 2, 1, 0., "");
      SLHAea_add(slha, "YD", 2, 2, sqrt2v*sminputs.mS, "s");
      SLHAea_add(slha, "YD", 2, 3, 0., "");
      SLHAea_add(slha, "YD", 3, 1, 0., "");
      SLHAea_add(slha, "YD", 3, 2, 0., "");
      SLHAea_add(slha, "YD", 3, 3, sqrt2v*sminputs.mBmB, "b");
      
      // Block SMINPUTS
      SLHAea_add_block(slha, "SMINPUTS");
      SLHAea_add(slha, "SMINPUTS", 1, sminputs.alphainv, "# alpha_em^-1(MZ)^MSbar");
      SLHAea_add(slha, "SMINPUTS", 2, sminputs.GF, "# G_mu [GeV^-2]");
      SLHAea_add(slha, "SMINPUTS", 3, sminputs.alphaS, "# alpha_s(MZ)^MSbar");
      SLHAea_add(slha, "SMINPUTS", 4, sminputs.mZ, "# m_Z(pole)");
      SLHAea_add(slha, "SMINPUTS", 5, sminputs.mBmB, "# m_b(m_b), MSbar");
      SLHAea_add(slha, "SMINPUTS", 6, sminputs.mT, "# m_t(pole)");
      SLHAea_add(slha, "SMINPUTS", 7, sminputs.mTau, "# m_tau(pole)");
      SLHAea_add(slha, "SMINPUTS", 8, sminputs.mNu3, "# m_nu_3");
      SLHAea_add(slha, "SMINPUTS", 11, sminputs.mE, "# m_e(pole)");
      SLHAea_add(slha, "SMINPUTS", 12, sminputs.mNu1, "# m_nu_1");
      SLHAea_add(slha, "SMINPUTS", 13, sminputs.mMu, "# m_muon(pole)");
      SLHAea_add(slha, "SMINPUTS", 14, sminputs.mNu2, "# m_nu_2");
      SLHAea_add(slha, "SMINPUTS", 21, sminputs.mD, "# m_d(2 GeV), MSbar");
      SLHAea_add(slha, "SMINPUTS", 22, sminputs.mU, "# m_u(2 GeV), MSbar");
      SLHAea_add(slha, "SMINPUTS", 23, sminputs.mS, "# m_s(2 GeV), MSbar");
      SLHAea_add(slha, "SMINPUTS", 24, sminputs.mCmC, "# m_c(m_c), MSbar");
      
      // And the W for good measure
      SLHAea_add(slha, "MASS", 24, sminputs.mW);
      
      // TODO: Chris Chang. Use the calculations of gggY
      double gggY1;
      double gggY2;
      std::tie(gggY1,gggY2) = getggY(sminputs.mT, *myPipe::Param["mY"], sminputs.alphaS, vev, *myPipe::Param["cY"]);
      SLHAea_add(slha, "DMINT", 3, gggY1);
      SLHAea_add(slha, "DMINT", 4, gggY2);

      // Retrieve any mass cuts
      static const Spectrum::mc_info mass_cut = myPipe::runOptions->getValueOrDef<Spectrum::mc_info>(Spectrum::mc_info(), "mass_cut");
      static const Spectrum::mr_info mass_ratio_cut = myPipe::runOptions->getValueOrDef<Spectrum::mr_info>(Spectrum::mr_info(), "mass_ratio_cut");
      
      // Construct the Spectrum object from the SLHAea inputs
      result = spectrum_from_SLHAea<Gambit::Models::MDMSMSimpleSpec, SLHAstruct>(slha,slha,mass_cut,mass_ratio_cut);
    }
    
    
    // Declaration: print spectrum out
    void fill_map_from_MDMSM_spectrum(std::map<std::string,double>&, const Spectrum&);
    
    void get_MDMSM_spectrum_as_map(std::map<std::string,double>& specmap)
    {
      namespace myPipe = Pipes::get_MDMSM_spectrum_as_map;
      const Spectrum& spec(*myPipe::Dep::MDMSM_spectrum);
      fill_map_from_MDMSM_spectrum(specmap, spec);
    }
    
    void fill_map_from_MDMSM_spectrum(std::map<std::string, double>& specmap, const Spectrum& spec)
    {
      /// Use SpectrumContents routines to automate
      static const SpectrumContents::MDMSM contents;
      static const std::vector<SpectrumParameter> required_parameters = contents.all_parameters();
      
      for(std::vector<SpectrumParameter>::const_iterator it = required_parameters.begin(); it != required_parameters.end(); ++it)
      {
        const Par::Tags        tag   = it->tag();
        const std::string      name  = it->name();
        const std::vector<int> shape = it->shape();
        
        // Scalar case
        if(shape.size()==1 and shape[0]==1)
        {
          std::ostringstream label;
          label << name <<" "<< Par::toString.at(tag);
          specmap[label.str()] = spec.get_HE().get(tag,name);
        }
        // Vector case
        else if(shape.size()==1 and shape[0]>1)
        {
          for(int i = 1; i<=shape[0]; ++i)
          {
            std::ostringstream label;
            label << name <<"_"<<i<<" "<< Par::toString.at(tag);
            specmap[label.str()] = spec.get_HE().get(tag,name,i);
          }
        }
        // Matrix case
        else if(shape.size()==2)
        {
          for(int i = 1; i<=shape[0]; ++i)
          {
            for(int j = 1; j<=shape[0]; ++j)
            {
              std::ostringstream label;
              label << name <<"_("<<i<<","<<j<<") "<<Par::toString.at(tag);
              specmap[label.str()] = spec.get_HE().get(tag,name,i,j);
            }
          }
        }
        // Deal with all other cases
        else
        {
          // ERROR
          std::ostringstream errmsg;
          errmsg << "Invalid parameter received while converting MDMSM_spectrum to map of strings!";
          errmsg << "Problematic parameter was: "<< tag <<", " << name << ", shape="<< shape;
          utils_error().forced_throw(LOCAL_INFO,errmsg.str());
        }
      }
    }
    
  }
  
}
