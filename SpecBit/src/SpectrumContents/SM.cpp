//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  Class defining the parameters that Spectrum
///  objects providing Standard Model parameters
///  (minus the Higgs sector) must provide.
///  Includes the SMInputs info plus pole masses
///
///  *********************************************
///
///  Authors:
///  <!-- add name and date if you modify -->
///
///  \author Ben Farmer
///          (benjamin.farmer@fysik.su.se)
///  \date 2016 Feb, 2019 Oct
///
///  \author Tomas Gonzalo
///          (tomas.gonzalo@monash.edu)
///  \date 2019 Oct
///
///  *********************************************

#include "gambit/SpecBit/RegisteredSpectra.hpp"

namespace Gambit {

  /// Only have to define the constructor
  SpectrumContents::SM::SM()
    : Contents("SM")
  {
     std::vector<int> scalar; // Empty vector means no indices   // i.e. get(Par::Tag, "name")
     std::vector<int> v3     = initVector(3);   // i.e. get(Par::Tag, "name", i)

     // TODO: I am not sure we need these
     //addParameter(Par::mass1, "u"   , v3, "MSBARMASS", 1); // TODO: Completely made up. Coordinate with SM_slha?
     //addParameter(Par::mass1, "d"   , v3, "MSBARMASS", 4);
     //addParameter(Par::mass1, "e-"  , v3, "MSBARMASS", 7);

     //addParameter(Par::mass1, "gamma", scalar, "MSBARMASS", 10);
     //addParameter(Par::mass1, "g"    , scalar, "MSBARMASS", 11);

     addParameter(Par::dimensionless, "alphainv", scalar, "SMINPUTS", 1); // TODO: Was 'alpha'. Needs consideration.
     addParameter(Par::imass2, "GF", scalar, "SMINPUTS", 2); // So far the only example of inverse mass dimension usage. Should be useful for other EFT couplings though.
     addParameter(Par::dimensionless, "alphaS",   scalar, "SMINPUTS", 3);

     addParameter(Par::Pole_Mass, "gamma", scalar, "MASS");
     addParameter(Par::Pole_Mass, "g"    , scalar, "MASS");
     addParameter(Par::Pole_Mass, "Z0"   , scalar, "SMINPUTS", 4);
     addParameter(Par::Pole_Mass, "W+"   , scalar, "MASS");

     addParameter(Par::Pole_Mass, "u_3" , scalar, "SMINPUTS", 6);
     addParameter(Par::mass1, "d_3" , scalar, "SMINPUTS", 5); // Not a pole mass but MSbar

     addParameter(Par::Pole_Mass, "nu_1", scalar, "SMINPUTS", 12);
     addParameter(Par::Pole_Mass, "nu_2", scalar, "SMINPUTS", 14);
     addParameter(Par::Pole_Mass, "nu_3", scalar, "SMINPUTS", 8);

     // SLHA SMINPUTS doesn't provide all those running masses. Only the light quarks:
     addParameter(Par::mass1, "u_1", scalar, "SMINPUTS", 22); // u
     addParameter(Par::mass1, "d_1", scalar, "SMINPUTS", 21); // d
     addParameter(Par::mass1, "d_2", scalar, "SMINPUTS", 23); // s
     addParameter(Par::mass1, "u_2", scalar, "SMINPUTS", 24); // c

     addParameter(Par::Pole_Mass, "e-_1"  , scalar, "SMINPUTS", 11);
     addParameter(Par::Pole_Mass, "e-_2" , scalar, "SMINPUTS", 13);
     addParameter(Par::Pole_Mass, "e-_3", scalar, "SMINPUTS", 7);


     addParameter(Par::Pole_Mixing, "sinW2", scalar, "SINTHETAW", 2);

     // Additional pole masses for SMInputs running masses
     addParameter(Par::Pole_Mass, "u_1", scalar, "MASS"); // u
     addParameter(Par::Pole_Mass, "d_1", scalar, "MASS"); // d
     addParameter(Par::Pole_Mass, "d_2", scalar, "MASS"); // s
     addParameter(Par::Pole_Mass, "u_2", scalar, "MASS"); // c
     addParameter(Par::Pole_Mass, "d_3", scalar, "MASS"); // b

     // Nearest flavour 'aliases' for the SM mass eigenstates
     addParameter(Par::Pole_Mass, "u", scalar, "MASS"); // u
     addParameter(Par::Pole_Mass, "d", scalar, "MASS"); // d
     addParameter(Par::Pole_Mass, "s", scalar, "MASS"); // s
     addParameter(Par::Pole_Mass, "c", scalar, "MASS"); // c
 
     addParameter(Par::Pole_Mass, "t" , scalar, "MASS");
     addParameter(Par::Pole_Mass, "b" , scalar, "MASS");

     addParameter(Par::Pole_Mass, "e-"  , scalar, "MASS");
     addParameter(Par::Pole_Mass, "mu-" , scalar, "MASS");
     addParameter(Par::Pole_Mass, "tau-", scalar, "MASS");

     // TODO: Considering adding CKM and PMNS here

  }

}
