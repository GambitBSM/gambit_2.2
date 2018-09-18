//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  Frontend source for the DirectDM backend.
///
///  *********************************************
///
///  Authors (add name and date if you modify):
///
///  \author Sanjay Bloor
///          (sanjay.bloor12@imperial.ac.uk)
///  \date 2018 Sep
///
///  *********************************************

#include "gambit/Backends/frontend_macros.hpp"
#include "gambit/Backends/frontends/DirectDM_2_0_1.hpp"

BE_NAMESPACE
{
  /* Convenience functions */

  /// Initialise a dictionary of Wilson Coefficients from a vector of string-double pairs,
  /// where the string is the name of the WC in DirectDM, and the double is the value at
  /// the input scale.
  void initialise_WC_dict(std::vector<std::pair<std::string, double>>>)
  {
    // e.g. in Python: dict1 = {'C61u' : 1./10000., 'C62u' : 1./10000., 'C61d' : 1./10000.}
  }

  /// Need to expose the following classes -- where DM_type can be equal to:
  /// "D": Dirac fermion, "M": Majorana fermion, "C": complex scalar, "R": real scalar DM. Default behaviour is DM_type="D".
  // class WC_3f(coeff_dict, DM_type)
  // class WC_4f(coeff_dict, DM_type)
  // class WC_5f(coeff_dict, DM_type)
  // class WC_EW(coeff_dict, Ychi, dchi, DM_type="D") -- only allows DM_type="D". Ychi, dchi are DM hypercharge and weak isospin.

  /// These classes all have the following methods:
  // run(mu_Lambda, muz=MZ): perform the one-loop RGE from scale mu_Lambda to muz (default MZ) and output a dictionary with the resulting Wilson coefficients. (Scales in GeV)

  /// match(DM_mass, mu_Lambda, RUN_EW=True): return a dictionary after performing the running and the matching to the theory with broken EW symmetry and five active quark flavors at scale MZ.
  // Setting the optional argument RUN_EW=False will switch off the RG evolution above the weak scale (default is RUN_EW=True).

  /// cNR(DM_mass, qvec, mu_Lambda, RGE=True, NLO=False, RUN_EW=True): return a dictionary containing the coefficients of the nuclear operators, c_i^N.
  // The two mandatory arguments are the DM mass (GeV) and the spatial momentum transfer qvec (GeV).
  // RGE=False will switch off the QCD and QED running (the default is RGE=True).
  // NLO=True will add the coherently enhanced NLO terms for the tensor operators (see Bishara et al. for the details); the default is NLO=False.
  // RUN_EW=False will switch off the RG evolution above the weak scale (default is RUN_EW=True).

}
END_BE_NAMESPACE
