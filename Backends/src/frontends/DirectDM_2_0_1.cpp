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
///  \date 2018 Sep, Oct
///
///  *********************************************

#include "gambit/Backends/frontend_macros.hpp"
#include "gambit/Backends/frontends/DirectDM_2_0_1.hpp"

BE_INI_FUNCTION
{
  // Empty ini function.
}
END_BE_INI_FUNCTION

/// TODO -
/// Set nuisance parameters for scans within directdm (BE_INI_FUNCTION)

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
    return relativistic_WCs;
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

    // Import Python class WC_EW from module DirectDM
    pybind11::object WC_EW = DirectDM.attr("WC_EW")(relativistic_WCs, Ychi, dchi, DM_type);
    
    // Obtain a dictionary of non-relativistic WCs, given the DM mass and the scale the Lagrangian is specified at.
    pybind11::dict cNRs = WC_EW.attr("_my_cNR")(mDM, scale);

    // Translate from Python dict to vec_strdbl_pairs structure
    pybind11::list cnrlistk = cNRs.attr("keys")();
    pybind11::list cnrlistv = cNRs.attr("values")();

    pybind11::str key;
    pybind11::float_ value;
    
    vec_strdbl_pairs nonrel_WCs;

    // Go through each key/value pair and append to the vector of non-rel WCs.
    for (unsigned int i=0; i<len(cnrlistk); i++)
    {
      key = cnrlistk[i];
      value = cnrlistv[i];
      nonrel_WCs.push_back(std::make_pair(std::string(key), (double)value));
    }

    return nonrel_WCs;
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

    // Python dictionary of non-relativistic Wilson Coefficients
    pybind11::dict cNRs;

    // Initialise Python class according to number of quark flavours specified from module DirectDM,
    // then obtain a dictionary of non-relativistic WCs, given the DM mass.
    if (scheme == 5)
    {
      pybind11::object WC_5f = DirectDM.attr("WC_5f")(relativistic_WCs, DM_type);
      cNRs = WC_5f.attr("_my_cNR")(mDM);
    }
    else if (scheme == 4)
    {
      pybind11::object WC_4f = DirectDM.attr("WC_4f")(relativistic_WCs, DM_type);
      cNRs = WC_4f.attr("_my_cNR")(mDM);
    }
    else if (scheme == 3)
    {
      pybind11::object WC_3f = DirectDM.attr("WC_3f")(relativistic_WCs, DM_type);
      cNRs = WC_3f.attr("_my_cNR")(mDM);
    }

    // Translate from Python dict to vec_strdbl_pairs structure
    pybind11::list cnrlistk = cNRs.attr("keys")();
    pybind11::list cnrlistv = cNRs.attr("values")();

    pybind11::str key;
    pybind11::float_ value;

    // Go through each key/value pair and append to the vector of non-rel WCs.
    for (unsigned int i=0; i<len(cnrlistk); i++)
    {
      key = cnrlistk[i];
      value = cnrlistv[i];
      nonrel_WCs.push_back(std::make_pair(std::string(key), (double)value));
    }

    return nonrel_WCs;
  }

}
END_BE_NAMESPACE
