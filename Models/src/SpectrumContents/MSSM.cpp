//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  Class defining the parameters that SubSpectrum
///  objects providing MSSM spectrum data must provide
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

#ifndef __mssmcontents_hpp__
#define __mssmcontents_hpp__

#include "gambit/Models/SpectrumContents/RegisteredSpectra.hpp"

namespace Gambit {

  /// Only have to define the constructor
  SpectrumContents::MSSM::MSSM()
  {
     setName("MSSM");

     addAllFrom(SM()); 

     // shape prototypes
     std::vector<int> scalar; // Empty vector, i.e. no indices, i.e.. get(Par::Tag, "name")
     std::vector<int> v2     = initVector(2);   // i.e. get(Par::Tag, "name", i)
     std::vector<int> v3     = initVector(3);   // "
     std::vector<int> v4     = initVector(4);   // "
     std::vector<int> v6     = initVector(6);   // "
     std::vector<int> m2x2   = initVector(2,2); // i.e. get(Par::Tag, "name", i, j)
     std::vector<int> m3x3   = initVector(3,3); // "
     std::vector<int> m4x4   = initVector(4,4); // "
     std::vector<int> m6x6   = initVector(6,6); // "

     //           tag,        name,   shape
     addParameter(Par::mass2, "BMu" , scalar, "BMu", 1); // TODO: Made this up
     addParameter(Par::mass2, "mHd2", scalar, "MSOFT", 1); // TODO: check order here, I forget which of mHu/d is mH1/2
     addParameter(Par::mass2, "mHu2", scalar, "MSOFT", 2);

     addParameter(Par::mass2, "mq2", m3x3, "MQ2"); // TODO: Pretty sure none of these are SLHA, so again have made up some blocks for them.
     addParameter(Par::mass2, "ml2", m3x3, "ML2");
     addParameter(Par::mass2, "md2", m3x3, "MD2");
     addParameter(Par::mass2, "mu2", m3x3, "MU2");
     addParameter(Par::mass2, "me2", m3x3, "ME2");

     addParameter(Par::mass1, "M1", scalar, "MSOFT", 1);
     addParameter(Par::mass1, "M2", scalar, "MSOFT", 2);
     addParameter(Par::mass1, "M3", scalar, "MSOFT", 3);
     addParameter(Par::mass1, "Mu", scalar, "HMIX", 1); // Is this \mu?
     addParameter(Par::mass1, "vu", scalar, "VEV", 1); // TODO: More made up stuff
     addParameter(Par::mass1, "vd", scalar, "VEV", 2);

     addParameter(Par::mass1, "TYd", m3x3, "TD"); // TODO: Peter check this. I think TD,TE etc are something to do with the trilinears in SLHA2, more like AD,AE etc, not this presumably Yukawa-related thing. I guess the SpecBit paper explains it. 
     addParameter(Par::mass1, "TYe", m3x3, "TE");
     addParameter(Par::mass1, "TYu", m3x3, "TU");
     addParameter(Par::mass1, "ad" , m3x3, "AD");
     addParameter(Par::mass1, "ae" , m3x3, "AE");
     addParameter(Par::mass1, "au" , m3x3, "AU");

     // EXTRAS! Kind of logical to always include these, without forcing users to calculate them themselves
     addParameter(Par::dimensionless, "tanbeta", scalar, "HMIX", 2); // DRBAR tanbeta
     //addParameter(Par::dimensionless, "tanbeta(mZ)", scalar); // i.e. the SLHA MINPAR value of tanbeta(mZ). Not yet a strict requirement, but highly recommended for wrappers to add it via override setters. TODO: Make this compulsory?
     addParameter(Par::mass2, "mA2" , scalar, "HMIX", 4);
     //

     addParameter(Par::dimensionless, "g1", scalar, "GAUGE", 1); // TODO: I forget if our g1,g2,g3 definitions are g',g,g3 like in SLHA. Check this. 
     addParameter(Par::dimensionless, "g2", scalar, "GAUGE", 2);
     addParameter(Par::dimensionless, "g3", scalar, "GAUGE", 3);

     addParameter(Par::dimensionless, "sinW2", scalar, "SINTHETAW", 1); // TODO: No SLHA definition for this, just make up somewhere for it

     addParameter(Par::dimensionless, "Yd", m3x3, "YU");
     addParameter(Par::dimensionless, "Yu", m3x3, "YD");
     addParameter(Par::dimensionless, "Ye", m3x3, "YE");

     addParameter(Par::Pole_Mass, "~g", scalar, "MASS");
     addParameter(Par::Pole_Mass, "~d",     v6, "MASS");
     addParameter(Par::Pole_Mass, "~u",     v6, "MASS");
     addParameter(Par::Pole_Mass, "~e-",    v6, "MASS");
     addParameter(Par::Pole_Mass, "~nu",    v3, "MASS");
     addParameter(Par::Pole_Mass, "~chi+",  v2, "MASS");
     addParameter(Par::Pole_Mass, "~chi0",  v4, "MASS");
     addParameter(Par::Pole_Mass, "h0",     v2, "MASS");
     addParameter(Par::Pole_Mass, "A0", scalar, "MASS");
     addParameter(Par::Pole_Mass, "H+", scalar, "MASS");
     addParameter(Par::Pole_Mass, "W+", scalar, "MASS");

     addParameter(Par::Pole_Mixing, "~d",    m6x6, "DSQMIX");
     addParameter(Par::Pole_Mixing, "~u",    m6x6, "USQMIX");
     addParameter(Par::Pole_Mixing, "~e-",   m6x6, "SELMIX");
     addParameter(Par::Pole_Mixing, "~nu",   m3x3, "SNUMIX");
     addParameter(Par::Pole_Mixing, "~chi0", m4x4, "NMIX");
     addParameter(Par::Pole_Mixing, "~chi-", m2x2, "UMIX"); // TODO: Does this naming really make the most sense?
     addParameter(Par::Pole_Mixing, "~chi+", m2x2, "VMIX"); 
     addParameter(Par::Pole_Mixing, "h0",    m2x2, "SCALARMIX"); // Non-SLHA: going with the FlexibleSUSY naming. Use these to compute ALPHA and other such SLHA variables when writing SLHA files.
     addParameter(Par::Pole_Mixing, "A0",    m2x2, "PSEUDOSCALARMIX");
     addParameter(Par::Pole_Mixing, "H+",    m2x2, "CHARGEMIX");
  }

}
#endif
