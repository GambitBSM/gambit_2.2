//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  Container for basic WIMP
///  particle properties
///
///  *********************************************
///
///  Authors (add name and date if you modify):
///
///  \author Ben Farmer
///          (b.farmer@imperial.ac.uk)
///  \date 2019 Sep
///
///  *********************************************

#ifndef __wimp_props_hpp__
#define __wimp_props_hpp__

#include <string>

namespace Gambit
{
    // Basic properties of generic WIMP
    struct WIMPprops
    {
      double mass;
      unsigned int spinx2;
      bool sc; // Self-conjugate?
      std::string name; // Name in the particle database
    };
}

#endif //__wimp_props_hpp__
