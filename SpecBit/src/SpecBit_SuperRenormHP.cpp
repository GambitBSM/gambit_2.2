//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  Functions of module SpecBit
///
///  SpecBit module functions related to the
///  super-renormalizable Higgs portal model.
///
///  *********************************************
///
///  Authors (add name and date if you modify):
///
///  \author Inigo Saez Casares
///          (inigo.saez_casares@ens-paris-saclay.fr)
///    \date 2020 March
///
///  *********************************************

#include <string>
#include <sstream>

#include "gambit/Elements/gambit_module_headers.hpp"

#include "gambit/Elements/spectrum.hpp"
#include "gambit/Utils/stream_overloads.hpp"
#include "gambit/Utils/util_macros.hpp"

#include "gambit/SpecBit/SpecBit_rollcall.hpp"
#include "gambit/SpecBit/SpecBit_helpers.hpp"
#include "gambit/SpecBit/QedQcdWrapper.hpp"
#include "gambit/Models/SimpleSpectra/SMHiggsSimpleSpec.hpp"
#include "gambit/Models/SimpleSpectra/SuperRenormHPSimpleSpec.hpp"
#include "gambit/SpecBit/model_files_and_boxes.hpp"

/* #include "gambit/SpecBit/ScalarSingletDM_Z2Spec.hpp" */
/* #include "gambit/SpecBit/SuperRenormHPSpec.hpp" */

// Switch for debug mode
//#define SPECBIT_DEBUG

namespace Gambit
{

  namespace SpecBit
  {
    using namespace LogTags;

    /// Get a (simple) Spectrum object wrapper for the SuperRenormHP model
    void get_SuperRenormHP_spectrum(Spectrum& result)
    {
      namespace myPipe = Pipes::get_SuperRenormHP_spectrum;
      const SMInputs& sminputs = *myPipe::Dep::SMINPUTS;

      // Initialise an object to carry the Singlet plus Higgs sector information
      Models::SuperRenormHPModel singletmodel;

      // quantities needed to fill container spectrum, intermediate calculations
      double alpha_em = 1.0 / sminputs.alphainv;
      double C = alpha_em * pi / (sminputs.GF * pow(2,0.5));
      double sinW2 = 0.5 - pow( 0.25 - C/pow(sminputs.mZ,2) , 0.5);
      double cosW2 = 0.5 + pow( 0.25 - C/pow(sminputs.mZ,2) , 0.5);
      double e = pow( 4*pi*( alpha_em ),0.5) ;

      // Higgs sector
      double mh = *myPipe::Param.at("mH");
      singletmodel.HiggsPoleMass = mh;

      double vev = 1. / sqrt(sqrt(2.)*sminputs.GF);
      singletmodel.HiggsVEV = vev;

      // Scalar singlet sector
      singletmodel.ScalarPoleMass = *myPipe::Param.at("mS");
      singletmodel.ScalarTheta = *myPipe::Param.at("theta");

      // Standard model
      singletmodel.sinW2 = sinW2;

      // gauge couplings
      singletmodel.g1 = sqrt(5/3) * e / sqrt(cosW2);
      singletmodel.g2 = e / sqrt(sinW2);
      singletmodel.g3   = pow( 4*pi*( sminputs.alphaS ),0.5) ;

      // Yukawas
      double sqrt2v = pow(2.0,0.5)/vev;
      singletmodel.Yu[0] = sqrt2v * sminputs.mU;
      singletmodel.Yu[1] = sqrt2v * sminputs.mCmC;
      singletmodel.Yu[2] = sqrt2v * sminputs.mT;
      singletmodel.Ye[0] = sqrt2v * sminputs.mE;
      singletmodel.Ye[1] = sqrt2v * sminputs.mMu;
      singletmodel.Ye[2] = sqrt2v * sminputs.mTau;
      singletmodel.Yd[0] = sqrt2v * sminputs.mD;
      singletmodel.Yd[1] = sqrt2v * sminputs.mS;
      singletmodel.Yd[2] = sqrt2v * sminputs.mBmB;

      // Create a SubSpectrum object to wrap the EW sector information
      Models::SuperRenormHPSimpleSpec singletspec(singletmodel);

      // Retrieve any mass cuts
      static const Spectrum::mc_info mass_cut = myPipe::runOptions->getValueOrDef<Spectrum::mc_info>(Spectrum::mc_info(), "mass_cut");
      static const Spectrum::mr_info mass_ratio_cut = myPipe::runOptions->getValueOrDef<Spectrum::mr_info>(Spectrum::mr_info(), "mass_ratio_cut");

      // We don't supply a LE subspectrum here; an SMSimpleSpec will therefore be automatically created from 'sminputs'
      result = Spectrum(singletspec,sminputs,&myPipe::Param,mass_cut,mass_ratio_cut);

    }

  } // end namespace SpecBit
} // end namespace Gambit

