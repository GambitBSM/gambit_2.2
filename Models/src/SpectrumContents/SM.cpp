//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  Class defining the parameters that SubSpectrum
///  objects providing Standard Model parameters
///  (minus the Higgs sector) must provide.
///
///  *********************************************
///
///  Authors:
///  <!-- add name and date if you modify -->
///
///  \author Ben Farmer
///          (benjamin.farmer@fysik.su.se)
///  \date 2016 Feb
///
///  *********************************************

#include "gambit/Models/SpectrumContents/RegisteredSpectra.hpp"

namespace Gambit {

  /// Only have to define the constructor
  SpectrumContents::SM::SM()
    : Contents("SM")
  {
     std::vector<int> scalar; // Empty vector means no indices   // i.e. get(Par::Tag, "name")
     std::vector<int> v3     = initVector(3);   // i.e. get(Par::Tag, "name", i)

     addParameter(Par::mass1, "u"   , v3, "MSBARMASS", 1); // TODO: Completely made up. Coordinate with SM_slha?
     addParameter(Par::mass1, "d"   , v3, "MSBARMASS", 4);
     addParameter(Par::mass1, "e-"  , v3, "MSBARMASS", 7);

     addParameter(Par::mass1, "gamma", scalar, "MSBARMASS", 10);
     addParameter(Par::mass1, "g"    , scalar, "MSBARMASS", 11);

     addParameter(Par::dimensionless, "alphainv", scalar, "SMINPUTS", 1); // TODO: Was 'alpha'. Needs consideration.
     addParameter(Par::dimensionless, "alphaS",   scalar, "SMINPUTS", 3);

     addParameter(Par::Pole_Mass, "gamma", scalar, "MASS");
     addParameter(Par::Pole_Mass, "g"    , scalar, "MASS");
     addParameter(Par::Pole_Mass, "Z0"   , scalar, "MASS");
     addParameter(Par::Pole_Mass, "W+"   , scalar, "MASS");

     addParameter(Par::Pole_Mass, "u_3" , scalar, "SMINPUTS", 6);
     addParameter(Par::Pole_Mass, "d_3" , scalar, "SMINPUTS", 5);

     addParameter(Par::Pole_Mass, "e-", v3, "MASS");
     addParameter(Par::Pole_Mass, "nu_1", scalar, "SMINPUTS", 12);
     addParameter(Par::Pole_Mass, "nu_2", scalar, "SMINPUTS", 14);
     addParameter(Par::Pole_Mass, "nu_3", scalar, "SMINPUTS", 8);


     addParameter(Par::Pole_Mixing, "sinW2", scalar, "SINTHETAW", 2);

     // Nearest flavour 'aliases' for the SM mass eigenstates
     addParameter(Par::Pole_Mass, "t" , scalar, "MASS");
     addParameter(Par::Pole_Mass, "b" , scalar, "MASS");

     addParameter(Par::Pole_Mass, "e-"  , scalar, "SMINPUTS", 11);
     addParameter(Par::Pole_Mass, "mu-" , scalar, "SMINPUTS", 13);
     addParameter(Par::Pole_Mass, "tau-", scalar, "SMINPUTS", 7);
  }

}
