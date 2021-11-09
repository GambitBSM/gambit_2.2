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

    void set_input_params(bool verbose, int NE_pd, int NT_pd, double eps)
    {
      py::module::import("acropolis.pprint").attr("verbose") = py::bool_(verbose);
      py::module::import("acropolis.cascade").attr("NE_pd") =  py::int_(NE_pd); // default: 150, fast: 75, aggresive: 30
      py::module::import("acropolis.cascade").attr("eps") = py::float_(eps); // default: 1e-3, fast: 1e-2, aggresive: 1e-1
      py::module::import("acropolis.nucl").attr("NT_pd") = py::int_(NT_pd); // default: 50, fast: 25, aggresive: 10
      py::module::import("acropolis.nucl").attr("eps") = py::float_(eps); // default: 1e-3, fast: 1e-2, aggresive: 1e-1
    }

    void abundance_photodissociation_decay(double* abundances_pre, double* abundances_post, double mass, double tau, double eta, double BR_el, double BR_ph)
    {
      #ifdef ACROPOLIS_DEBUG
        std::cout << "[ACROPOLIS] Invoking 'DecayModel' with (mass, tau, T0, eta, BR_el, BR_ph) = ";
        std::cout << "( " << mass << " " << tau << " " << T0 << " " << eta << " " << BR_el << " " << BR_ph << " )" << std::endl;
      #endif

      // Initialise the model
      py::object mod = AC_models.attr("DecayModel")(mass, tau, T0, eta, BR_el, BR_ph);

      // Get the initial abundances of the isotopes n, p, H2, H3, He3, He4, Li6, Li7, Be7
      pyArray_dbl initial_abundances(9, abundances_pre);

      // Reshape the numpy array [shape (9,)] into a 1D-matrix [shape (9,1)]
      initial_abundances = initial_abundances.attr("reshape")(9,1);
 
      // Replace the internal initial abundance matrix of the 'InputInterface' with the content of 'intial_abundances'
      mod.attr("_sII").attr("set_bbn_abundances")(initial_abundances);
 
      // Run the disintegration and compute the final abundances
      pyArray_dbl final_abundances = mod.attr("run_disintegration")();

      // Write the results into 'abundances_post'
      for (int i=0; i != final_abundances.size(); ++i)
      {
        *(abundances_post+i) = *(final_abundances.data()+i);
      }
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
