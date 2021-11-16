//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  Frontend for Acropolis 1.2.1 backend
///
///  *********************************************
///
///  Authors (add name and date if you modify):
///
///  \author Patrick Stoecker
///          (patrick.stoecker@kit.edu)
///  \date 2021 Oct
///
///  \author Tomas Gonzalo
///          (gonzalo@physik.rwth-aachen.de)
///  \date 2021 Oct
///
///  *********************************************

#include "gambit/Backends/frontend_macros.hpp"
#include "gambit/Backends/frontends/Acropolis_1_2_1.hpp"

#include "gambit/Utils/numerical_constants.hpp"

#include <algorithm>

//#define ACROPOLIS_DEBUG

#ifdef HAVE_PYBIND11

  // Convenience functions (definitions)
  BE_NAMESPACE
  {
    namespace py = pybind11;

    using pyArray_dbl = py::array_t<double>;

    // Static variables to store important modules and submodules
    static py::module AC;
    static py::module AC_models;

    // Reference temperature [in MeV] at which the initial number density of the decaying particle is defined.
    static double T0; 

    // Translate the results of AlterBBN (ratioH) into the abundances (Y0) that are assumed by ACROPOLIS.
    //
    // This function is adapted from tools/create_sm_abundance_file.c  on the ACROPOLIS git repository.
    // https://github.com/hep-mh/acropolis.git [added at commit 59d76e9]
    //
    // (The only change is the implicit shift in the index: AlterBBN starts counting at 1, other people would start at 0)
    void ratioH_to_Y0(double ratioH[], double Y0[], int NY)
    {
      // Handle the special case 'p'
      Y0[1] = ratioH[1];
      // Handle the secial case 'He4'
      Y0[5] = ratioH[5]/4;

      // All other isotopes are normalised with respect to 'p'
      for ( int i = 0; i < NY; i++ )
      {
        // Skip the special cases from above
        if ( i == 1 || i == 5 ) continue;
        Y0[i] = ratioH[i]*Y0[1];
      }

      // Revert the decays of 'H3' and 'Be7'
      Y0[4] = Y0[4] - ratioH[3]*Y0[1]; // H3
      Y0[7] = Y0[7] - ratioH[8]*Y0[1]; // Be7

      // Perform the neutron decay
      // TODO: This happens after photodissociation, not now
      //Y0[1] = Y0[1] + Y0[0];           // p
      //Y0[0] = 0.;                      // n
    }

    // Reverse translation
    void Y0_to_ratioH(double ratioH[], double Y0[], int NY)
    {
      // Handle the special case 'p'
      ratioH[1] = Y0[1];
      // Handle the secial case 'He4'
      ratioH[5] = Y0[5]*4;

      // All other isotopes are normalised with respect to 'p'
      for ( int i = 0; i < NY; i++ )
      {
        // Skip the special cases from above
        if ( i == 1 || i == 5 ) continue;
        ratioH[i] = Y0[i]/Y0[1];
      }

      // Perform the 'post-BBN' decays of 'H3' and 'Be7'
      // (c.f. plots/plot_annih_pwave_1MeV_ee.py of the ACROPOLIS git repository)
      // TODO: This decays are already computed by ACROPOLIS
      //ratioH[4] = ratioH[4] + ratioH[3]; // H3 -> He3
      //ratioH[7] = ratioH[7] + ratioH[8]; // Be7 -> li7

    }

    void set_input_params(bool verbose, int NE_pd, int NT_pd, double eps)
    {
      py::module::import("acropolis.pprint").attr("verbose") = py::bool_(verbose);
      py::module::import("acropolis.cascade").attr("NE_pd") =  py::int_(NE_pd); // default: 150, fast: 75, aggresive: 30
      py::module::import("acropolis.cascade").attr("eps") = py::float_(eps); // default: 1e-3, fast: 1e-2, aggresive: 1e-1
      py::module::import("acropolis.nucl").attr("NT_pd") = py::int_(NT_pd); // default: 50, fast: 25, aggresive: 10
      py::module::import("acropolis.nucl").attr("eps") = py::float_(eps); // default: 1e-3, fast: 1e-2, aggresive: 1e-1
    }

    void abundance_photodissociation_decay(double* ratioH_pre, double* cov_ratioH_pre, double* ratioH_post, double* cov_ratioH_post, double mass, double tau, double eta, double BR_el, double BR_ph, int niso)
    {
      #ifdef ACROPOLIS_DEBUG
        std::cout << "[ACROPOLIS] Invoking 'DecayModel' with (mass, tau, T0, eta, BR_el, BR_ph) = ";
        std::cout << "( " << mass << " " << tau << " " << T0 << " " << eta << " " << BR_el << " " << BR_ph << " )" << std::endl;
      #endif

      // Initialise the model
      py::object mod = AC_models.attr("DecayModel")(mass, tau, T0, eta, BR_el, BR_ph);

      // Get the initial abundances (in the 'AlterBBN basis') of the isotopes n, p, H2, H3, He3, He4, Li6, Li7, Be7
      pyArray_dbl Y0_pre(niso, ratioH_pre);

      // Translate the results of AlterBBN ('ratioH') into pure abundances ('Y0') needed by ACROPOLIS.
      ratioH_to_Y0(ratioH_pre, Y0_pre.mutable_data(), niso);

      // Reshape the numpy array [shape (niso,)] into a 1D-matrix [shape (niso,1)]
      Y0_pre = Y0_pre.attr("reshape")(niso,1);

      // Replace the internal initial abundance matrix of the 'InputInterface' with the content of 'intial_abundances'
      mod.attr("_sII").attr("set_bbn_abundances")(Y0_pre);

      // Run the disintegration and compute the final abundances
      py::dict result = mod.attr("run_disintegration")();

      // Get the transfer matrix
      pyArray_dbl transfer_matrix = py::cast<pyArray_dbl>(result["transfer_matrix"]);

      // Pure abundances post acropolis
      double Y0_post[niso];

      // Write the results into 'ratioH_post' and 'cov_ratioH_post'
      for (int i=0; i != niso; ++i)
      {
        Y0_post[i] = 0.0;
        for (int j=0; j < niso; ++j)
        {
          *(Y0_post+i) += *(transfer_matrix.data()+i*niso+j) * *(Y0_pre.data()+j);
          // TODO: Set the covariances to 0 for now. Might expand later
          *(cov_ratioH_post+i*niso+j) = 0.0;
          //for (int k=0; k < niso; ++k) for (int l=0; l < niso; ++l)
          //{
          //  *(cov_ratioH_post+i*niso+j) += *(transfer_matrix.data()+i*niso+k) *  *(cov_ratioH_pre+k*niso+l) * *(transfer_matrix.data()+j*niso+l);
          //}
        }
      }

      // Transalte the results from pure abundances to ratioH
      Y0_to_ratioH(ratioH_post, Y0_post, niso);
    }

  }
  END_BE_NAMESPACE

#endif

// Initialisation function (definition)
BE_INI_FUNCTION
{

  #ifdef HAVE_PYBIND11

    static bool first_point = true;

    // Enter this scope only for the first point
    if (first_point)
    {
      first_point = false;

      // Save the submodule "models" (for later)
      AC = Gambit::Backends::Acropolis_1_2_1::Acropolis;
      std::string module_name = AC.attr("__name__").cast<std::string>();
      AC_models = py::module::import( (module_name + ".models").c_str() );

      // Set the reference temperature to 10 MeV
      T0 = 10.0;
    }
  
  #endif
}
END_BE_INI_FUNCTION
