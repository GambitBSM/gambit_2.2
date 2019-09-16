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


    /*! \brief Two-dimensional integration container for bilinear interpolation and bicubic splines.
     */

    // XrayInterpolator2D class: Provides a 2-D interpolation container based on the gsl library.
    // Can be declared static for efficiency & easy one-time initialisation of interpolating functions.
    class XrayInterpolator2D
    {
      public:
        // Overloaded class creators for the XrayInterpolator class using the init function below.
        XrayInterpolator2D();
        XrayInterpolator2D(std::string file, std::string type);
        XrayInterpolator2D(std::string file);
        XrayInterpolator2D& operator=(XrayInterpolator2D&&);
        // Destructor.
        ~XrayInterpolator2D();
        // Delete copy constructor and assignment operator to avoid shallow copies.
        XrayInterpolator2D(const XrayInterpolator2D&) = delete;
        XrayInterpolator2D operator=(const XrayInterpolator2D&) = delete;
        // Routine to access interpolated values.
        double interpolate(double x, double y);
        // Routine to check if a point is inside the interpolating box.
        bool is_inside_box(double x, double y);
      private:
        // Initialiser for the XrayInterpolator2D class.
        void init(std::string file, std::string type);
        // The gsl objects for the interpolating functions that need to be available to the class routines.
        gsl_interp_accel *x_acc;
        gsl_interp_accel *y_acc;
        gsl_spline2d *spline;
        // Upper and lower "x" and "y" values available to the interpolating function.
        double x_lo, y_lo, x_up, y_up;
    };

    // Move assignment operator
    XrayInterpolator2D& XrayInterpolator2D::operator=(XrayInterpolator2D&& interp)
    {
      if(this != &interp)
      {
        std::swap(x_acc,interp.x_acc);
        std::swap(y_acc,interp.y_acc);
        std::swap(spline,interp.spline);
        std::swap(x_lo,interp.x_lo);
        std::swap(x_up,interp.x_up);
        std::swap(y_lo,interp.y_lo);
        std::swap(y_up,interp.y_up);
      }
      return *this;
    }

    // Destructor
    XrayInterpolator2D::~XrayInterpolator2D()
    {
      gsl_spline2d_free (spline);
      gsl_interp_accel_free (x_acc);
      gsl_interp_accel_free (y_acc);
    }

    // Initialiser for the XrayInterpolator class.
    void XrayInterpolator2D::init(std::string file, std::string type)
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
      tab.setcolnames("x", "y", "z");
      // Initialise gsl interpolation routine.
      // Get unique entries of "x" and "y" for the grid and grid size.
      std::vector<double> x_vec = tab["x"];
      sort(x_vec.begin(), x_vec.end());
      x_vec.erase(unique(x_vec.begin(), x_vec.end()), x_vec.end());
      int nx = x_vec.size();
      std::vector<double> y_vec = tab["y"];
      sort(y_vec.begin(), y_vec.end());
      y_vec.erase(unique(y_vec.begin(), y_vec.end()), y_vec.end());
      int ny = y_vec.size();
      int n_grid_pts = tab["z"].size();

      if (nx*ny != n_grid_pts)
      {
        DarkBit_error().raise(LOCAL_INFO, "ERROR! The number of grid points ("+std::to_string(n_grid_pts)+") for XrayInterpolator2D does not equal the number of unique 'x' and 'y' values ("+std::to_string(nx)+" and "+std::to_string(ny)+")!\n       Check formatting of the file: '"+file+"'.");
      };

      const double* x = &x_vec[0];
      const double* y = &y_vec[0];
      // Allocate memory for "z" values array in gsl format
      double* z = (double*) malloc(nx * ny * sizeof(double));

      if (type == "bicubic")
      {
        spline = gsl_spline2d_alloc(gsl_interp2d_bicubic, nx, ny);
      }
      else if (type == "bilinear")
      {
        spline = gsl_spline2d_alloc(gsl_interp2d_bilinear, nx, ny);
      }
      else
      {
        DarkBit_error().raise(LOCAL_INFO, "ERROR! Interpolation type '"+type+"' not known to class XrayInterpolator2D.\n       Available types: 'bilinear' and 'bicubic'.");
      };

      x_acc = gsl_interp_accel_alloc();
      y_acc = gsl_interp_accel_alloc();

      // Determine first and last "x" and "y" values and grid step size.
      x_lo = x_vec.front();
      x_up = x_vec.back();
      y_lo = y_vec.front();
      y_up = y_vec.back();
      double x_delta = (x_up-x_lo) / (nx-1);
      double y_delta = (y_up-y_lo) / (ny-1);

      // Intialise grid.
      for (int i = 0; i < n_grid_pts; i++)
      {
        // Determine appropriate indices for the grid points.
        double temp = (tab["x"][i]-x_lo) / x_delta;
        int ind_x = (int) (temp+0.5);
        temp = (tab["y"][i]-y_lo) / y_delta;
        int ind_y = (int) (temp+0.5);

        //std::cout << ind_x << "/" << nx-1 << " " << tab["x"][i] << " vs " << x[ind_x] << " " << ind_y << "/" << ny-1 << " " << tab["y"][i] << " vs " << y[ind_y] << std::endl;

        gsl_spline2d_set(spline, z, ind_x, ind_y, tab["z"][i]);
      };
        gsl_spline2d_init (spline, x, y, z, nx, ny);
    };

    // Overloaded class creators for the Axion}Interpolator class using the init function above.
    XrayInterpolator2D::XrayInterpolator2D(std::string file, std::string type) { init(file, type); };
    XrayInterpolator2D::XrayInterpolator2D(std::string file) { init(file, "bilinear"); };
    XrayInterpolator2D::XrayInterpolator2D() {};

    // Routine to access interpolated values.
    double XrayInterpolator2D::interpolate(double x, double y) { return gsl_spline2d_eval(spline, x, y, x_acc, y_acc); };

    // Routine to check if a point is inside the interpolating box.
    bool XrayInterpolator2D::is_inside_box(double x, double y) { return ((x >= x_lo) && (x <= x_up) && (y >= y_lo) && (y <= y_up)); };


    ///////////////////////////////////////
    //         Xray likelihoods          //
    ///////////////////////////////////////

    void compute_lnL_Xray_toy(double& result)
    {
      using namespace Pipes::compute_lnL_Xray_toy;

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

      double logm = log10(*Param["mass"]) + 9; // In "DecayingDM_photon", the mass is given in GeV. Need to convert it into eV
      double tau = *Param["lifetime"];   // lifetime is already in untis of s. No tranformation needed.
      double frac = *Param["fraction"];

      if (logm <= xlim.first || logm >= xlim.second)
      {
        // Bound can only be applied if log10(mass) is within the range of the table.
        result = 0.0;
      }
      else
      {
        double tau_bound = pow(10.,WISPy_bound.interpolate(logm));
        bool within_bound = ((1./frac)*exp(t_universe/tau)*tau < tau_bound);
        result = (within_bound ? -9.0 : 0.0);
      }
    }
  }
}
