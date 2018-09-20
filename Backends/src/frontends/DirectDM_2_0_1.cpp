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

#include <pybind11/embed.h>

BE_INI_FUNCTION
{
  // Empty ini function.
}
END_BE_INI_FUNCTION

BE_NAMESPACE
{
  /* Convenience functions */

  /// Initialise a dictionary of relativistic Wilson Coefficients from a vector of
  /// string-double pairs, where the string is the name of the WC in DirectDM, and
  /// the double is the value at the input scale.
  pybind11::dict initialise_WC_dict(vec_strdbl_pairs& wilsonCoeffs)
  {
    pybind11::dict relativistic_WCs;
    // e.g. in Python: dict1 = {'C61u' : 1./10000., 'C62u' : 1./10000., 'C61d' : 1./10000.}

    // Loop through all Wilson Coefficients and assign them to the Python dict.
    for (auto it = wilsonCoeffs.begin(); it != wilsonCoeffs.end(); ++it)
    {
      relativistic_WCs[pybind11::cast(it->first)] = pybind11::cast(it->second);
    }
    pybind11::print(relativistic_WCs);
    return relativistic_WCs;
  }

  /// Get Wilson Coefficients at 2 GeV in a quark flavour matching scheme.
  /// Requires a dictionary of relatavistic WCs, the DM mass, an integer specifying the number
  /// of quark flavours to match onto, i.e. the 3, 4 or 5 quark flavour scheme, and
  /// the DM type -- "D" for Dirac fermion; "M" for Majorana fermion; "C" for complex scalar; "R" for real scalar.
  vec_strdbl_pairs get_NR_WCs_flav(pybind11::dict& relativistic_WCs, double& mDM, int& scheme, std::string& DM_type)
  {
    // Only load up 3, 4, 5 flavour scheme.
    if (scheme != 3 || scheme != 4 || scheme != 5)
    {
      backend_error().raise(LOCAL_INFO, "DirectDM quark flavour matching scheme must be for "
        "3, 4 or 5 quark flavors.");
    }
    vec_strdbl_pairs nonrel_WCs;

    // Import the Python class "WC_3f" from "directdm"
    /// Todo - replace by directdm object loaded up by BE system (loaded_python_backends)
    pybind11::object wc_3f = pybind11::module::import("directdm").attr("WC_3f");
    pybind11::object wc_4f = pybind11::module::import("directdm").attr("WC_4f");
    pybind11::object wc_5f = pybind11::module::import("directdm").attr("WC_5f");

    // Initialise DirectDM object according to number of flavours specified
    if (scheme == 5)
    {
      pybind11::object WCs = wc_5f(relativistic_WCs, DM_type);
    }


    return nonrel_WCs;
  }

  /// Get Wilson Coefficients at 2 GeV from the SM unbroken phase.
  /// Requires a dictionary of relatavistic WCs, , the DM mass, dchi is the dimension of the DM SU2 representation,
  /// Ychi is the DM hypercharge such that Q = I^3 + Y/2, scale is the scale the Lagrangian is defined at, and
  /// the DM type -- "D" for Dirac fermion; "M" for Majorana fermion; "C" for complex scalar; "R" for real scalar.
  vec_strdbl_pairs get_NR_WCs_EW(pybind11::dict& relativistic_WCs,  double& mDM, double& dchi, double& Ychi, double& scale, std::string& DM_type)
  {
    // S.B. 19/09: currently only Dirac supported
    if (DM_type != "D")
    {
      backend_error().raise(LOCAL_INFO, "DirectDM at unbroken scale currenly only supports Dirac DM.");
    }
    vec_strdbl_pairs nonrel_WCs;

    // Import Python class WC_EW from module directdm
    pybind11::object ddm = pybind11::module::import("directdm");
    pybind11::object WC_EW = ddm.attr("WC_EW")(relativistic_WCs, Ychi, dchi, DM_type);
    pybind11::object cNRs = WC_EW.attr("_my_cNR")(mDM, scale);

    pybind11::print(cNRs);


    return nonrel_WCs;
  }

  /// Need to expose the following classes -- where DM_type can be equal to:
  /// "D": Dirac fermion, "M": Majorana fermion, "C": complex scalar, "R": real scalar DM. Default behaviour is DM_type="D".
  // class WC_3f(coeff_dict, DM_type)
  // class WC_4f(coeff_dict, DM_type)
  // class WC_5f(coeff_dict, DM_type)
  // class WC_EW(coeff_dict, Ychi, dchi, DM_type="D") -- only allows DM_type="D". Ychi, dchi are DM hypercharge and dimension of the SU(2) representation.

  /// These classes all have the following methods:
  // run(mu_Lambda, muz=MZ): perform the one-loop RGE from scale mu_Lambda to muz (default MZ) and output a dictionary with the resulting Wilson coefficients. (Scales in GeV)

  /// match(DM_mass, mu_Lambda, RUN_EW=True): return a dictionary after performing the running and the matching to the theory with broken EW symmetry and five active quark flavors at scale MZ.
  // Setting the optional argument RUN_EW=False will switch off the RG evolution above the weak scale (default is RUN_EW=True).

  /// cNR(DM_mass, qvec, mu_Lambda, RGE=True, NLO=False, RUN_EW=True): return a dictionary containing the coefficients of the nuclear operators, c_i^N.
  /// _my_cNR(DM_mass, mu_Lambda, RGE=True, NLO=False, RUN_EW=True): return a dictionary containing the coefficients of the nuclear operators, c_i^N.
  // The two mandatory arguments are the DM mass (GeV) and the spatial momentum transfer qvec (GeV).
  // RGE=False will switch off the QCD and QED running (the default is RGE=True).
  // NLO=True will add the coherently enhanced NLO terms for the tensor operators (see Bishara et al. for the details); the default is NLO=False.
  // RUN_EW=False will switch off the RG evolution above the weak scale (default is RUN_EW=True).

}
END_BE_NAMESPACE
