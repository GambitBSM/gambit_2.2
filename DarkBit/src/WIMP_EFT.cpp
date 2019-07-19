//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  Module functions associated with creating
///  and translating WIMP-nucleon and WIMP-quark 
///  effective operator couplings from GAMBIT
///  ModelParameters. Functions which compute
///  these EFT couplings for specific "UV" models 
///  live in DarkBit sources files named after those
///  models.
///
///  Includes module functions to compute
///  non-relativistic operator couplings from
///  relativistic ones using DirectDM.
///
///  *********************************************
///
///  Authors (add name and date if you modify):
///
///  \author Ben Farmer
///          (b.farmer@imperial.ac.uk)
///  \date 2019 Jul
///
///  *********************************************

#include "gambit/Elements/gambit_module_headers.hpp"
#include "gambit/DarkBit/DarkBit_rollcall.hpp"
#include "gambit/DarkBit/DarkBit_utils.hpp"

//#define DARKBIT_DEBUG

namespace Gambit
{

  namespace DarkBit
  {
 
    //////////////////////////////////////////////////////////////////////////
    //
    //   Translation of NREO ModelParameters into NREO_DM_nucleon_couplings
    //
    //////////////////////////////////////////////////////////////////////////

    void NREO_couplings_from_parameters(NREO_DM_nucleon_couplings& NREO_couplings)
    {
       using namespace Pipes::NREO_couplings_from_parameters;
       NREO_couplings = NREO_DM_nucleon_couplings(Param); // Constructor takes care of the parameter copying for us
    }

    //////////////////////////////////////////////////////////////////////////
    //
    //   Translation of DD_couplings into NREO_DM_nucleon_couplings
    //
    //////////////////////////////////////////////////////////////////////////

    void NREO_from_DD_couplings(NREO_DM_nucleon_couplings& NREO_couplings)
    {
       using namespace Pipes::NREO_from_DD_couplings;
       DM_nucleon_couplings ddc = *Dep::DD_couplings;

       // TODO! I have not been able to find the exact conventions
       // used in DDcalc vs the NREO model. I think it is just this:
       // c0 = 0.5*(cp+cn)
       // c1 = 0.5*(cp-cn)
       // so that 
       // cp = c0 + c1
       // cn = c0 - c1
       // Change if needed!
    
       // Compute non-zero isospin basis couplings from DM_nucleon_couplings entries
       // TODO: I also did this from memory, should check I got the operator numbers right
       NREO_couplings.c0[1] = 0.5*(ddc.gps + ddc.gns);
       NREO_couplings.c1[1] = 0.5*(ddc.gps - ddc.gns);
       NREO_couplings.c0[4] = 0.5*(ddc.gpa + ddc.gna);
       NREO_couplings.c1[4] = 0.5*(ddc.gpa - ddc.gna);
    }

    /* Non-relativistic Wilson Coefficients, model independent */

    /// Obtain the non-relativistic Wilson Coefficients from a set of model
    /// specific relativistic Wilson Coefficients from DirectDM in the flavour
    /// matching scheme (default 5 flavours). NR WCs defined at 2 GeV.
    void DD_nonrel_WCs_flavscheme(NREO_DM_nucleon_couplings &result)
    {
      using namespace Pipes::DD_nonrel_WCs_flavscheme;

      // Number of quark flavours used for matching (default 5)
      int scheme = runOptions->getValueOrDef<int>(5,"flavs");

      // Obtain spin of DM particle, plus identify whether DM is self-conjugate
      double mDM = *Dep::mwimp;
      unsigned int sDM  = *Dep::spinwimpx2;
      bool is_SC = *Dep::wimp_sc;

      // Set DM_type based on the spin and & conjugacy of DM
      std::string DM_type;

      // Fermion case
      if (sDM == 1) { is_SC ? DM_type == "M" : DM_type == "D"; }
      // Scalar
      else if (sDM == 0) { is_SC ? DM_type == "R" : DM_type == "C"; }

      // Relativistic Wilson Coefficients
      map_str_dbl relativistic_WCs = *Dep::DD_rel_WCs;

      // Get non-relativistic coefficients
      result = BEreq::get_NR_WCs_flav(relativistic_WCs, mDM, scheme, DM_type);
    }

    /// Obtain the non-relativistic Wilson Coefficients from a set of model
    /// specific relativistic Wilson Coefficients from DirectDM in the
    /// unbroken SM phase. NR WCs defined at 2 GeV.
    void DD_nonrel_WCs_EW(NREO_DM_nucleon_couplings &result)
    {
      using namespace Pipes::DD_nonrel_WCs_EW;

      // Specify the scale that the Lagrangian is defined at
      double scale = runOptions->getValue<double>("scale");
      // Hypercharge of DM
      double Ychi = runOptions->getValue<double>("Ychi");
      // SU(2) dimension of DM
      double dchi = runOptions->getValue<int>("dchi");

      // Obtain spin of DM particle, plus identify whether DM is self-conjugate
      double mDM = *Dep::mwimp;
      unsigned int sDM  = *Dep::spinwimpx2;
      bool is_SC = *Dep::wimp_sc;

      // Set DM_type based on the spin and & conjugacy of DM
      std::string DM_type;

      // Set DM_type based on the spin and & conjugacy of DM
      // Fermion case: set DM_type to Majorana or Dirac
      if (sDM == 1) { is_SC ? DM_type = "M" : DM_type = "D"; }
      // Scalar case: set DM type to real or complex
      else if (sDM == 0) { is_SC ? DM_type = "R" : DM_type = "C"; }

      // Relativistic Wilson Coefficients
      map_str_dbl relativistic_WCs = *Dep::DD_rel_WCs;

      // Get non-relativistic coefficients
      /// TODO - How to get hypercharge and SU(2) dimension for these fields!?
      /// Currently just comes from the YAML file. GUM? Process Catalogue?
      result = BEreq::get_NR_WCs_EW(relativistic_WCs, mDM, dchi, Ychi, scale, DM_type);
    }

  }
}
