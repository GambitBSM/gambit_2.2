//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  Class defining the parameters that SubSpectrum
///  objects providing the Standard Model Higgs
///  sector parameters must provide
///
///  *********************************************
///
///  Authors:
///  <!-- add name and date if you modify -->
///
///  \author Ben Farmer
///          (benjamin.farmer@imperial.ac.uk)
///  \date 2016 Feb, 2019 Oct
///
///  *********************************************

#include "gambit/SpecBit/SpectrumContents/RegisteredSpectra.hpp"

namespace Gambit {

  /// Only have to define the constructor
  SpectrumContents::SMHiggs::SMHiggs()
    : Contents("SMHiggs")
  {
     std::vector<int> scalar; // Empty vector, i.e. no indices, i.e.. get(Par::Tag, "name")
 
     addParameter(Par::mass1,     "vev",  scalar, "VEVS", 1);
     addParameter(Par::Pole_Mass, "h0_1", scalar, "MASS");
  }

}
