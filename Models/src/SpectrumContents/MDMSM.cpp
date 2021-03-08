//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  Class defining the parameters that SubSpectrum
///  objects providing MDMSM
///  spectrum data must provide.
///
///  Authors (add name and date if you modify):    
///       *** Automatically created by GUM ***     
///                                                
///  \author The GAMBIT Collaboration             
///  \date 04:53PM on March 08, 2021
///                                                
///  ********************************************* 

#ifndef __MDMSM_contents_hpp__
#define __MDMSM_contents_hpp__

#include "gambit/Models/SpectrumContents/RegisteredSpectra.hpp"

namespace Gambit
{
  SpectrumContents::MDMSM::MDMSM()
  {
    setName("MDMSM");
    
    std::vector<int> scalar = initVector(1); // i.e. get(Par::Tag, "name")
    std::vector<int> m3x3  = initVector(3,3); // i.e. get(Par::Tag, "name", i, j)
    
    addParameter(Par::dimensionless, "gchi", scalar, "DMINT", 1);
    addParameter(Par::dimensionless, "cY", scalar, "DMINT", 2);
    addParameter(Par::dimensionless, "gggY1", scalar, "DMINT", 3);
    addParameter(Par::dimensionless, "gggY2", scalar, "DMINT", 4);
    addParameter(Par::mass1, "vev", scalar, "VEVS", 1);
    addParameter(Par::dimensionless, "g1", scalar, "GAUGE", 1);
    addParameter(Par::dimensionless, "g2", scalar, "GAUGE", 2);
    addParameter(Par::dimensionless, "g3", scalar, "GAUGE", 3);
    addParameter(Par::dimensionless, "sinW2", scalar, "SINTHETAW", 1);
    addParameter(Par::dimensionless, "Yd", m3x3, "YD", 1);
    addParameter(Par::dimensionless, "Yu", m3x3, "YU", 1);
    addParameter(Par::dimensionless, "Ye", m3x3, "YE", 1);
    addParameter(Par::Pole_Mass, "h0_1", scalar, "MASS", 25);
    addParameter(Par::Pole_Mass, "~chi", scalar, "MASS", 52);
    addParameter(Par::Pole_Mass, "Y", scalar, "MASS", 99902);
    
  } // namespace Models
} // namespace Gambit
#endif
