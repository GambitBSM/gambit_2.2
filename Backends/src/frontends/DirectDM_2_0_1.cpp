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
///  \author Ben Farmer
///          (b.farmer@imperial.ac.uk)
///  \date 2019 Jul
///
///  *********************************************

#include "gambit/Backends/frontend_macros.hpp"
#include "gambit/Backends/frontends/DirectDM_2_0_1.hpp"
#include "gambit/Backends/backend_types/DDCalc.hpp" // Using DDCalc (frontend) container objects

#include <pybind11/stl.h>

BE_INI_FUNCTION
{
  // Empty ini function.
}
END_BE_INI_FUNCTION

/// TODO -
/// Set nuisance parameters for scans within directdm (BE_INI_FUNCTION)

BE_NAMESPACE
{
  /* Simple local helper function to copy map of coefficients into GAMBIT (DDCalc) container */
  NREO_DM_nucleon_couplings copy_couplings_to_NREO_container(const map_str_dbl& nonrel_WCs)
  {
    // Copy coefficients into GAMBIT container object (converting to isoscalar/isovector basis in the process)
    // Assuming conventions of 1203.3542, i.e.
    // c0 = cp + cn
    // c1 = cp - cn
    NREO_DM_nucleon_couplings NRWCs;

    for(int OpCoeff=1; OpCoeff<=15; OpCoeff++)
    {
       std::stringstream sp;
       std::stringstream sn;
       sp<<"cNR"<<OpCoeff<<"p";
       sp<<"cNR"<<OpCoeff<<"n";
       auto icp = nonrel_WCs.find(sp.str());
       auto icn = nonrel_WCs.find(sn.str());
       double cp=0;
       double cn=0;
       if(icp!=nonrel_WCs.end()) cp=(icp->second);
       if(icn!=nonrel_WCs.end()) cn=(icn->second);
       NRWCs.c0[OpCoeff] = cp + cn;
       NRWCs.c1[OpCoeff] = cp - cn;
    } 
    return NRWCs;
  }

  /* Convenience functions */

  /// Get Wilson Coefficients at 2 GeV from the SM unbroken phase.
  /// Requires a dictionary of relatavistic WCs, , the DM mass, dchi is the dimension of the DM SU2 representation,
  /// Ychi is the DM hypercharge such that Q = I^3 + Y/2, scale is the scale the Lagrangian is defined at, and
  /// the DM type -- "D" for Dirac fermion; "M" for Majorana fermion; "C" for complex scalar; "R" for real scalar.
  NREO_DM_nucleon_couplings get_NR_WCs_EW(map_str_dbl& relativistic_WCs,  double& mDM, double& dchi, double& Ychi, double& scale, std::string& DM_type)
  {
    // S.B. 19/09/18: currently only Dirac supported
    if (DM_type != "D")
    {
      backend_error().raise(LOCAL_INFO, "DirectDM at unbroken scale currenly only supports Dirac DM.");
    }

    // Import Python class WC_EW from module DirectDM
    pybind11::object WC_EW = DirectDM.attr("WC_EW")(relativistic_WCs, Ychi, dchi, DM_type);
    
    // Obtain a dictionary of non-relativistic WCs, given the DM mass and the scale the Lagrangian is specified at.
    pybind11::dict cNRs = WC_EW.attr("_my_cNR")(mDM, scale);
    
    // Cast python dictionary to C++ type known to GAMBIT    
    map_str_dbl nonrel_WCs = cNRs.cast<map_str_dbl>(); 

    // Copy coefficients into GAMBIT container object (converting to isoscalar/isovector basis in the process)
    return copy_couplings_to_NREO_container(nonrel_WCs);
  }

  /// Get Wilson Coefficients at 2 GeV in a quark flavour matching scheme.
  /// Requires a dictionary of relatavistic WCs, the DM mass, an integer specifying the number
  /// of quark flavours to match onto, i.e. the 3, 4 or 5 quark flavour scheme, and
  /// the DM type -- "D" for Dirac fermion; "M" for Majorana fermion; "C" for complex scalar; "R" for real scalar.
  NREO_DM_nucleon_couplings get_NR_WCs_flav(map_str_dbl& relativistic_WCs, double& mDM, int& scheme, std::string& DM_type)
  {
    // Only load up 3, 4, 5 flavour scheme.
    if (scheme != 3 || scheme != 4 || scheme != 5)
    {
      backend_error().raise(LOCAL_INFO, "DirectDM quark flavour matching scheme must be for "
        "3, 4 or 5 quark flavors.");
    }

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

    // Cast python dictionary to C++ type known to GAMBIT    
    map_str_dbl nonrel_WCs = cNRs.cast<map_str_dbl>(); 

    // Copy coefficients into GAMBIT container object (converting to isoscalar/isovector basis in the process)
    return copy_couplings_to_NREO_container(nonrel_WCs);
  }

}
END_BE_NAMESPACE
