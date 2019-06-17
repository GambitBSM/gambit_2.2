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

namespace Gambit
{

  /// Only have to define the constructor
  SpectrumContents::SM_slha::SM_slha()
  {
     setName("SM_slha");

     std::vector<int> scalar = initVector(1);   // i.e. get(Par::Tag, "name")
     std::vector<int> v3     = initVector(3);   // i.e. get(Par::Tag, "name", i)

     //addParameter(Par::mass1, "u"   , v3);
     //addParameter(Par::mass1, "d"   , v3);
     //addParameter(Par::mass1, "e-"  , v3);
     //
     // SLHA SMINPUTS doesn't provide all those running masses. Only the light quarks:
     addParameter(Par::mass1, "u_1", scalar, "SMINPUTS", 22); // u
     addParameter(Par::mass1, "d_1", scalar, "SMINPUTS", 21); // d
     addParameter(Par::mass1, "d_2", scalar, "SMINPUTS", 23); // s

     // doesn't provide these either, although they are just zero of course
     //addParameter(Par::mass1, "gamma", scalar);
     //addParameter(Par::mass1, "g"    , scalar);

     // And we don't provide these either since they cannot be provided at the same scale.
     //addParameter(Par::dimensionless, "alpha" , scalar);
     //addParameter(Par::dimensionless, "alphaS", scalar);

     addParameter(Par::Pole_Mass, "gamma", scalar, "MASS");
     addParameter(Par::Pole_Mass, "g"    , scalar, "MASS");
     addParameter(Par::Pole_Mass, "Z0"   , scalar, "SMINPUTS", 4);
     addParameter(Par::Pole_Mass, "W+"   , scalar, "MASS");

     addParameter(Par::Pole_Mass, "u_3" , scalar, "SMINPUTS", 6); // t
     addParameter(Par::Pole_Mass, "d_3" , scalar, "SMINPUTS", 5); // b (technically mb(mb)^{MSbar}, not pole)

     addParameter(Par::Pole_Mass, "e-", v3, "MASS");
     addParameter(Par::Pole_Mass, "nu_1", scalar, "SMINPUTS", 12);
     addParameter(Par::Pole_Mass, "nu_2", scalar, "SMINPUTS", 14);
     addParameter(Par::Pole_Mass, "nu_3", scalar, "SMINPUTS", 8);

     addParameter(Par::Pole_Mixing, "sinW2", scalar, "SINTHETAW", 2); // TODO: Need to sort out exact definition of this. 

     // Nearest flavour 'aliases' for the SM mass eigenstates
     addParameter(Par::Pole_Mass, "t" , scalar, "MASS");
     addParameter(Par::Pole_Mass, "b" , scalar, "MASS");

     addParameter(Par::Pole_Mass, "e-"  , scalar, "SMINPUTS", 11);
     addParameter(Par::Pole_Mass, "mu-" , scalar, "SMINPUTS", 13);
     addParameter(Par::Pole_Mass, "tau-", scalar, "SMINPUTS", 7);
  }

}
