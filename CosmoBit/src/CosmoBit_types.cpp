//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  Source code for types for module CosmoBit.
///  For instructions on adding new types, see
///  the corresponding header.
///
///  *********************************************
///
///  Authors (add name and date if you modify):
///
///  \author Patrick Stoecker
///          (stoecker@physik.rwth-aachen.de)
///  \date 2017 Nov
///
///  \author Selim Hotinli
///  \date 2018 Jan
///
///  \author Janina Renk
///          (janina.renk@fysik.su.se)
///  \date 2018 Oct
///
///  \author Sebastian Hoof
///          (hoof@uni-goettingen.de)
///  \date 2020 Mar
///
///  \author Pat Scott
///          (pat.scott@uq.edu.au)
///  \date 2020 Apr
///
///  *********************************************

#include <string>
#include <iostream>

#include "gambit/CosmoBit/CosmoBit_types.hpp"
#include "gambit/CosmoBit/CosmoBit_utils.hpp"
#include "gambit/Utils/numerical_constants.hpp"

namespace Gambit
{
  namespace CosmoBit
  {

    SM_time_evo::SM_time_evo(double t0, double tf, double N_t, double Neff_SM) : grid_size(N_t), t_grid(N_t), T_evo(N_t), Tnu_evo(N_t), H_evo(N_t), H_int(N_t)
    {

      // check if implemented routines are valid for given initial time
      if(t0 < 1e3)
      {
        std::ostringstream err;
        err << "Requested initial time for evolution of Temperature & Hubble rate for SM for Temperatures t_initial was "<< t0<<". Implemented routines are only valid after e+/- annihilation (t > 10^3). Aborting now.";
        CosmoBit_error().raise(LOCAL_INFO, err.str());
      }

      // initialize time grid in log space
      double Delta_logt = (log(tf) - log(t0))/(grid_size-1);
      for (int jj = 0; jj<grid_size; ++jj)
      {
        t_grid[jj] = exp(log(t0) + jj*Delta_logt);
      }
      double g_star_SM = 2.+2.*7./8.*Neff_SM*pow(4./11.,4./3.); // contribution from photons & neutrinos with Neff = 3.046

      // factor needed to calculate temperature evolution. For details see definition of functions set_T_evo(),.. in CosmoBit_types.hpp header
      factor_T_evo = 1./sqrt(2.*sqrt(8.*pi*pi*pi*GN_SI *pow(kB_SI,4.)*g_star_SM/90./c_SI/c_SI/pow(hP_SI/2./pi*c_SI,3.)))*kB_eV_over_K/1e3;
      factor_Tnu_evo = pow(Neff_SM/3.,1./4.)* pow(4./11.,1./3.)*factor_T_evo; // = T_nu/T_gamma * factor_T_evo
      factor_HT_evo = sqrt(8.*pi*GN_SI/3.*pi*pi/30.*g_star_SM/pow(hP_SI/2./pi*c_SI,3.))*(1e6*eV_to_J*eV_to_J)/c_SI;

      // set time evolution of photon T, neutrino T, and Hubble
      // rate based on above calculated evolution factors 
      // for time grid 't_grid'
      set_T_evo();
      set_Tnu_evo();
      set_Ht_evo();
    }

    void SM_time_evo::calc_H_int()
    {
      /*
      This calculates int(y(x') dx', {x', x0, x}) for each x in x_grid (with x0 = x_grid[0]), using a simple trapezoidal integration.
      This function explicitly assumes that the logarithms of the values in x_grid are equally spaced.
      This is very simplified and designed to work fast in case of the Hubble rate. Can go wrong with other functions, so it should
      really only be used in this context (or if you exactly know what you are doing..).
      */
      std::valarray<double> g_grid(grid_size);
      double Delta_z = (log(t_grid[grid_size-1]) - log(t_grid[0]))/(grid_size - 1);

      g_grid = t_grid*H_evo;

      H_int[0] = 0.0;
      H_int[1] = 0.5*Delta_z*(g_grid[0] + g_grid[1]);

      double cumsum = g_grid[1];
      for (int i = 2; i<grid_size; ++i)
      {
          H_int[i] = 0.5*Delta_z*(g_grid[0] + g_grid[i] + 2.0*cumsum);
          cumsum = cumsum + g_grid[i];
      }
    }

    // return Parametrised_ps members A_s, n_s, r, and N_pivot as str to double map for printing
    map_str_dbl Parametrised_ps::get_parametrised_ps_map()
    {
      map_str_dbl result;
      result["ln10A_s"] = ln10A_s;
      result["n_s"] = n_s;
      result["r"] = r;
      result["N_pivot"] = N_pivot;

      return result;
    };

    void Primordial_ps::fill_k(double *k_array, int len)
    {
      std::vector<double> K(k_array, k_array+len);
      k = std::move(K);
      vec_size = len;
    }

    void Primordial_ps::fill_P_s(double *P_s_array, int len)
    {
      std::vector<double> ps(P_s_array, P_s_array+len);
      P_s = std::move(ps);
    }

    void Primordial_ps::fill_P_s_iso(double *P_s_iso_array, int len)
    {
      std::vector<double> psi(P_s_iso_array, P_s_iso_array+len);
      P_s_iso = std::move(psi);
    }

    void Primordial_ps::fill_P_t(double *P_t_array, int len)
    {
      std::vector<double> pt(P_t_array, P_t_array+len);
      P_t = std::move(pt);
    }
  }
}
