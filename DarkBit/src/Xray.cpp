///  GAMBIT: Global and Modular BSM Inference Tool
///  *********************************************
///  \file
///
///  Xray likelihoods for DarkBit.
///
///  *********************************************
///
///  Authors (add name and date if you modify):
///
///  \author Patrick Stöcker
///          (stoecker@physik.rwth-aachen.de)
///  \date 2019 Sep
///
///  *********************************************

#include <algorithm>
#include <cmath>
#include <vector>
#include <string>
#include <iostream>
#include <sstream>

#include <gsl/gsl_math.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_sf_trig.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>
#include <gsl/gsl_histogram.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>

#include "gambit/Elements/gambit_module_headers.hpp"
#include "gambit/Utils/util_functions.hpp"
#include "gambit/Utils/ascii_table_reader.hpp"
#include "gambit/Utils/statistics.hpp"
#include "gambit/Utils/numerical_constants.hpp"
#include "gambit/DarkBit/DarkBit_rollcall.hpp"
#include "gambit/DarkBit/DarkBit_utils.hpp"

namespace Gambit
{
  namespace DarkBit
  {
    /////////////////////////////////////////////////////////////////
    //      Auxillary functions and classes for interpolation      //
    /////////////////////////////////////////////////////////////////

    /*! \brief Generic one-dimensional integration container for linear interpolation and cubic splines.
     */

    // XrayInterpolator class: Provides a general 1-D interpolation container based on the gsl library.
    // Can be declared static for efficiency & easy one-time initialisation of interpolating functions.
    // (This is the twin sibling of AxionInterpolator in DarkBit/src/Axions.cpp)
    class XrayInterpolator
    {
      public:
        // Overloaded class creators for the XrayInterpolator class using the init function below.
        XrayInterpolator();
        XrayInterpolator(const std::vector<double> x, const std::vector<double> y, std::string type);
        XrayInterpolator(const std::vector<double> x, const std::vector<double> y);
        XrayInterpolator(std::string file, std::string type);
        XrayInterpolator(std::string file);
        XrayInterpolator& operator=(XrayInterpolator&&);
        // Destructor.
        ~XrayInterpolator();
        // Delete copy constructor and assignment operator to avoid shallow copies.
        XrayInterpolator(const XrayInterpolator&) = delete;
        XrayInterpolator operator=(const XrayInterpolator&) = delete;
        // Routine to access interpolated values.
        double interpolate(double x);
        // Routine to access upper and lower boundaries of available data.
        double lower();
        double upper();
      private:
        // Initialiser for the XrayInterpolator class.
        void init(std::string file, std::string type);
        void init(const std::vector<double> x, const std::vector<double> y, std::string type);
        // The gsl objects for the interpolating functions.
        gsl_interp_accel *acc;
        gsl_spline *spline;
        // Upper and lower boundaries available for the interpolating function.
        double lo;
        double up;
    };

    // Default constructor.
    XrayInterpolator::XrayInterpolator() {};

    // Initialiser for the XrayInterpolator class.
    void XrayInterpolator::init(const std::vector<double> x, const std::vector<double> y, std::string type)
    {
      int pts = x.size();
      // Get first and last value of the "x" component.
      lo = x.front();
      up = x.back();
      acc = gsl_interp_accel_alloc ();
      if (type == "cspline")
      {
        spline = gsl_spline_alloc (gsl_interp_cspline, pts);
      }
      else if (type == "linear")
      {
        spline = gsl_spline_alloc (gsl_interp_linear, pts);
      }
      else
      {
        DarkBit_error().raise(LOCAL_INFO, "ERROR! Interpolation type '"+type+"' not known to class XrayInterpolator.\n       Available types: 'linear' and 'cspline'.");
      };

      gsl_spline_init (spline, &x[0], &y[0], pts);
    };

    // Overloaded class creators for the XrayInterpolator class using the init function above.
    XrayInterpolator::XrayInterpolator(const std::vector<double> x, const std::vector<double> y, std::string type) { init(x, y, type); };
    XrayInterpolator::XrayInterpolator(const std::vector<double> x, const std::vector<double> y) { init(x, y, "linear"); };

    // Initialiser for the XrayInterpolator class.
    void XrayInterpolator::init(std::string file, std::string type)
    {
      // Check if file exists.
      if (not(Utils::file_exists(file)))
      {
        DarkBit_error().raise(LOCAL_INFO, "ERROR! File '"+file+"' not found!");
      } else {
        logger() << LogTags::debug << "Reading data from file '"+file+"' and interpolating it with '"+type+"' method." << EOM;
      };
      // Read numerical values from data file.
      ASCIItableReader tab (file);
      tab.setcolnames("x", "y");

      // for (int idx=1; idx < tab["x"].size(); idx++)
      // {
      //   std::cout << "x[" << idx << "] = " << tab["x"][idx] << "; dx = " << tab["x"][idx] -tab["x"][idx-1] << std::endl;
      //   if (tab["x"][idx] -tab["x"][idx-1] <= 0.0) std::cout << "OH NO" << std::endl;
      // }

      init(tab["x"],tab["y"],type);
    };

    // Overloaded class creators for the XrayInterpolator class using the init function above.
    XrayInterpolator::XrayInterpolator(std::string file, std::string type) { init(file, type); };
    XrayInterpolator::XrayInterpolator(std::string file) { init(file, "linear"); };

    // Move assignment operator
    XrayInterpolator& XrayInterpolator::operator=(XrayInterpolator&& interp)
    {
      if(this != &interp)
      {
        std::swap(acc,interp.acc);
        std::swap(spline,interp.spline);
        std::swap(lo,interp.lo);
        std::swap(up,interp.up);
      }
      return *this;
    }

    // Destructor
    XrayInterpolator::~XrayInterpolator()
    {
      gsl_spline_free (spline);
      gsl_interp_accel_free (acc);
    }

    // Routine to access interpolated values.
    double XrayInterpolator::interpolate(double x) { return gsl_spline_eval(spline, x, acc); };

    // Routines to return upper and lower boundaries of interpolating function
    double XrayInterpolator::lower() { return lo; };
    double XrayInterpolator::upper() { return up; };

    ////////////////////////////////////////////////////
    //               Xray likelihoods                 //
    // -- based on likehoods for sterile neutrinos -- //
    ////////////////////////////////////////////////////

    void compute_lnL_Xray_WISPy(double& result)
    {
      using namespace Pipes::compute_lnL_Xray_WISPy;

      static XrayInterpolator WISPy_bound;
      static bool first = true;
      static std::pair<double,double> xlim;

      if (first)
      {
        WISPy_bound = std::move(XrayInterpolator(GAMBIT_DIR "/DarkBit/data/WISPy_bound.dat","linear"));
        xlim.first = WISPy_bound.lower();
        xlim.second = WISPy_bound.upper();
        first = false;
      }

      const double t_universe = 4.32e17; // Age of the universe in seconds (https://www.physicsoftheuniverse.com/numbers.html)

      double logm = log10(*Param["mass"]) + 9; // In "DecayingDM_general", the mass is given in GeV. Need to convert it into eV
      double tau = *Param["lifetime"];   // lifetime is already in untis of s. No tranformation needed.
      double frac = *Param["fraction"];
      double BR = *Param["BR"];

      if (logm <= xlim.first || logm >= xlim.second)
      {
        // Bound can only be applied if log10(mass) is within the range of the table.
        result = 0.0;
      }
      else
      {
        double tau_bound = pow(10.,WISPy_bound.interpolate(logm));
        bool excluded = ((1./frac)*exp(t_universe/tau)*BR*tau < tau_bound);
        result = (excluded ? -9.0 : 0.0);
      }
    }

    void compute_lnL_Xray_Integral(double& result)
    {
      using namespace Pipes::compute_lnL_Xray_Integral;

      static XrayInterpolator sin2_2t_bound;
      static bool first = true;
      static std::pair<double,double> xlim;

      if (first)
      {
        sin2_2t_bound = std::move(XrayInterpolator(GAMBIT_DIR "/DarkBit/data/Integral_sterile_nu_bound.dat","linear"));
        xlim.first = sin2_2t_bound.lower();
        xlim.second = sin2_2t_bound.upper();
        first = false;
      }

      const double t_universe = 4.32e17; // Age of the universe in seconds (https://www.physicsoftheuniverse.com/numbers.html)

      double mass = *Param["mass"] * 1e6; // In "DecayingDM_general", the mass is given in GeV. Need to convert it into keV
      double tau = *Param["lifetime"];   // lifetime is already in untis of s. No tranformation needed.
      double frac = *Param["fraction"];
      double BR = *Param["BR"];

      if (mass <= xlim.first || mass >= xlim.second)
      {
        // Bound can only be applied if the mass is within the range of the table.
        result = 0.0;
      }
      else
      {
        double sin2_2t = sin2_2t_bound.interpolate(mass);
        double tau_bound = 2./1.36038e-32 * pow(1e10*sin2_2t,-1.)*pow(mass,-5.);
        bool excluded = ((1./frac)*exp(t_universe/tau)*BR*tau < tau_bound);
        result = (excluded ? -9.0 : 0.0);
      }
    }

    void compute_lnL_Xray_M31(double& result)
    {
      using namespace Pipes::compute_lnL_Xray_M31;

      static XrayInterpolator sin2_2t_bound;
      static bool first = true;
      static std::pair<double,double> xlim;

      if (first)
      {
        sin2_2t_bound = std::move(XrayInterpolator(GAMBIT_DIR "/DarkBit/data/M31_sterile_nu_bound.dat","linear"));
        xlim.first = sin2_2t_bound.lower();
        xlim.second = sin2_2t_bound.upper();
        first = false;
      }

      const double t_universe = 4.32e17; // Age of the universe in seconds (https://www.physicsoftheuniverse.com/numbers.html)

      double mass = *Param["mass"] * 1e6; // In "DecayingDM_general", the mass is given in GeV. Need to convert it into keV
      double tau = *Param["lifetime"];   // lifetime is already in untis of s. No tranformation needed.
      double frac = *Param["fraction"];
      double BR = *Param["BR"];

      if (mass <= xlim.first || mass >= xlim.second)
      {
        // Bound can only be applied if the mass is within the range of the table.
        result = 0.0;
      }
      else
      {
        double sin2_2t = sin2_2t_bound.interpolate(mass);
        double tau_bound = 2./1.36038e-32 * pow(1e10*sin2_2t,-1.)*pow(mass,-5.);
        bool excluded = ((1./frac)*exp(t_universe/tau)*BR*tau < tau_bound);
        result = (excluded ? -9.0 : 0.0);
      }
    }

    void compute_lnL_Xray_NuSTAR(double& result)
    {
      using namespace Pipes::compute_lnL_Xray_NuSTAR;

      static XrayInterpolator sin2_2t_bound;
      static bool first = true;
      static std::pair<double,double> xlim;

      if (first)
      {
        sin2_2t_bound = std::move(XrayInterpolator(GAMBIT_DIR "/DarkBit/data/NuSTAR_sterile_nu_bound.dat","linear"));
        xlim.first = sin2_2t_bound.lower();
        xlim.second = sin2_2t_bound.upper();
        first = false;
      }

      const double t_universe = 4.32e17; // Age of the universe in seconds (https://www.physicsoftheuniverse.com/numbers.html)

      double mass = *Param["mass"] * 1e6; // In "DecayingDM_general", the mass is given in GeV. Need to convert it into keV
      double tau = *Param["lifetime"];   // lifetime is already in untis of s. No tranformation needed.
      double frac = *Param["fraction"];
      double BR = *Param["BR"];

      if (mass <= xlim.first || mass >= xlim.second)
      {
        // Bound can only be applied if the mass is within the range of the table.
        result = 0.0;
      }
      else
      {
        double sin2_2t = sin2_2t_bound.interpolate(mass);
        double tau_bound = 2./1.36038e-32 * pow(1e10*sin2_2t,-1.)*pow(mass,-5.);
        bool excluded = ((1./frac)*exp(t_universe/tau)*BR*tau < tau_bound);
        result = (excluded ? -9.0 : 0.0);
      }
    }
  }
}
