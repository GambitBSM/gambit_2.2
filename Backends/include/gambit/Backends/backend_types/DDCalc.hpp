//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  Helper types for DDCalc backend.
///
///  *********************************************
///
///  Authors (add name and date if you modify):
///
///  \author Pat Scott
///          p.scott@imperial.ac.uk
///  \date 2016 May
///
///  *************************

#ifndef __DDCalc_types_hpp__
#define __DDCalc_types_hpp__

namespace Gambit
{

  // Container for dark matter - nucleon couplings
  struct DM_nucleon_couplings
  {
    double gps;
    double gns;
    double gpa;
    double gna;
  };

  struct DD_coupling_container
  {
    int coeff_structure;  // Simple integer to tell DDCalc which effective operators to set the WIMP object up with.
    DM_nucleon_couplings          DM_nucleon_coeffs; // Corresponds to int = 1. Direct DM-nucleon interactions.
    std::map<std::string,double>  DD_nonrel_WCs;     // Corresponds to int = 2. Effective non-relativistic DM-quark Wilson coefficients.
  };

}

#endif /* defined __DDCalc_types_hpp__ */
