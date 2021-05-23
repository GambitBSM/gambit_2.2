//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  ColliderBit classes for holding a 
///  collection of 1D/2D interpolators
///
///  *********************************************
///
///  Authors (add name and date if you modify):
///
///  \author Anders Kvellestad
///          (anders.kvellestad@fys.uio.no)
///  \date 2021 May
///
///  *********************************************

#include <vector>
#include <map>
#include <string>

#include <gsl/gsl_spline.h>
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>

#include "gambit/Utils/ascii_table_reader.hpp"
#include "gambit/Utils/util_functions.hpp"


namespace Gambit
{

  namespace ColliderBit
  {

    /// A class for holding a collection of 1D interpolators, created from reading a tabulated ascii file. 
    /// - The first column is taken as the x value. 
    /// - For each remaining column a 1D interpolation function f(x) is created
    class interp1d_collection
    {
      public:

        // ==== Member variables ====

        string collection_name;
        string file_name;
        string x_name;
        vector<string> interpolator_names;

        vector<gsl_spline*> splines;
        vector<gsl_interp_accel*> x_accels;
        
        double x_min;
        double x_max;
        size_t n_interpolators;


        // ==== Methods ====

        // Constructor
        interp1d_collection(const string collection_name_in, const string filename_in, const vector<string> colnames_in)
        {
          // Check if file exists.
          if (not(Utils::file_exists(filename_in)))
          {
            ColliderBit_error().raise(LOCAL_INFO, "ERROR! File '" + filename_in + "' not found!");
          }

          // Read numerical values from data file.
          ASCIItableReader tab(filename_in);

          // Check that there's more than one column
          if (tab.getncol() < 2)
          {
            ColliderBit_error().raise(LOCAL_INFO, "ERROR! Less than two columns found in the input file '" + filename_in + "'."); 
          }

          // Check that the number of columns match the number of column names
          if (tab.getncol() != colnames.size())
          {
            ColliderBit_error().raise(LOCAL_INFO, "ERROR! Mismatch between number of columns and number of column names."); 
          }

          // Set the column names
          tab.setcolnames(colnames);

          // Save some names
          collection_name = collection_name_in;
          file_name = filename_in;
          x_name = colnames_in.at(0);
          interpolator_names = {colnames_in.begin() + 1, colnames_in.end()};

          // Number of interpolators
          n_interpolators = interpolator_names.size();

          // 
          // Create one GSL spline for each interpolator
          // 
          vector<double> x_data = tab[x_name];
          int nx = x_data.size();

          // Get first and last value of the "x" component.
          x_min = x_data.front();
          x_max = x_data.back();

          for (size_t interp_index = 0; interp_index < n_interpolators; ++interp_index)
          {
            interp_name = interpolator_names[interp_index];
            vector<double> interp_data = tab[interp_name];

            int n_points = interp_data.size();
            if (nx != n_points)
            {
              ColliderBit_error().raise(LOCAL_INFO, "ERROR! The number of data points ("+std::to_string(n_points)+") does not agree with the number of 'x' points ("+std::to_string(nx)+") for the interpolator '"+interp_name+"'.\n Check formatting of the file: '"+file+"'.");
            }

            // Initialize a gsl_spline pointer and a gsl_interp_accel pointer and store them
            gsl_interp_accel* acc = gsl_interp_accel_alloc();
            x_accels.push_back(acc);

            gsl_spline* sp = gsl_spline_alloc(gsl_interp_linear, pts);  // For now only using linear interpolation
            gsl_spline_init(sp, &x_data[0], &interp_data[0], pts);
            splines.push_back(sp);
          }
        }

        // Destructor
        ~interp1d_collection()
        {
          for (size_t interp_index = 0; interp_index < n_interpolators; ++interp_index)
          {
            gsl_spline_free(splines[interp_index]);
            gsl_interp_accel_free(x_accels[interp_index]);
          }
        }

        // Evaluate a given interpolation
        double eval(double x, size_t interp_index)
        {
          return gsl_spline_eval(splines[interp_index], x, x_accels[interp_index]);
        }
        
        // If there is only one interpolator, we can evaluate it without specifying the index
        double eval(double x)
        {
          if (n_interpolators != 1)
          {
              ColliderBit_error().raise(LOCAL_INFO, "ERROR! This interp1d_collection instance contains more than one interpolator, so the interpolator index must be specified.");
          }
          return eval(x, 0);
        }

        // Check if point is inside interpolation range
        bool is_inside_range(double x)
        {
          return ((x >= x_min) && (x <= x_max));
        }

    };


    /// A class for holding a collection of 2D interpolators, created from reading a tabulated ascii file. 
    /// - The first two columns are taken to be the x and y grid points. 
    /// - For each remaining column a 2D interpolation function f(x,y) is created
    /// - Note that GLS assumes the points are ordered according to increasing x (first) and y (second) values, 
    ///   i.e. an ordering like (x0,y0), (x0,y1), (x0,y2), ... (x1,y0), (x1,y1), (x1,y2), ...
    ///
    class interp2d_collection
    {
      public:

        // ==== Member variables ====

        string collection_name;
        string file_name;
        string x_name;
        string y_name;
        vector<string> interpolator_names;

        vector<gsl_spline2d*> splines;
        vector<gsl_interp_accel*> x_accels;
        vector<gsl_interp_accel*> y_accels;
        
        double x_min;
        double x_max;
        double y_min;
        double y_max;

        size_t n_interpolators;


        // ==== Methods ====

        // Constructor
        interp2d_collection(const string collection_name_in, const string filename_in, const vector<string> colnames_in)
        {
          // Check if file exists.
          if (not(Utils::file_exists(filename_in)))
          {
            ColliderBit_error().raise(LOCAL_INFO, "ERROR! File '" + filename_in + "' not found!");
          }

          // Read numerical values from data file.
          ASCIItableReader tab(filename_in);

          // Check that there's more than two columns
          if (tab.getncol() < 3)
          {
            ColliderBit_error().raise(LOCAL_INFO, "ERROR! Less than three columns found in the input file '" + filename_in + "'."); 
          }

          // Check that the number of columns match the number of column names
          if (tab.getncol() != colnames.size())
          {
            ColliderBit_error().raise(LOCAL_INFO, "ERROR! Mismatch between number of columns and number of column names."); 
          }

          // Set the column names
          tab.setcolnames(colnames);

          // Save some names
          collection_name = collection_name_in;
          file_name = filename_in;
          x_name = colnames_in.at(0);
          y_name = colnames_in.at(1);
          interpolator_names = {colnames_in.begin() + 2, colnames_in.end()};

          // Number of interpolators
          n_interpolators = interpolator_names.size();

          // 
          // Create one GSL spline for each interpolator
          // 

          // Get unique entries of "x" and "y" for the grid and grid size.
          std::vector<double> x_vec = tab[x_name];
          sort(x_vec.begin(), x_vec.end());
          x_vec.erase(unique(x_vec.begin(), x_vec.end()), x_vec.end());
          int nx = x_vec.size();
          x_min = x_vec.front();
          x_max = x_vec.back();

          std::vector<double> y_vec = tab[y_name];
          sort(y_vec.begin(), y_vec.end());
          y_vec.erase(unique(y_vec.begin(), y_vec.end()), y_vec.end());
          int ny = y_vec.size();
          y_min = y_vec.front();
          y_max = y_vec.back();

          for (size_t interp_index = 0; interp_index < n_interpolators; ++interp_index)
          {
            interp_name = interpolator_names[interp_index];
            vector<double> interp_data = tab[interp_name];

            int n_points = interp_data.size();

            if (nx * ny != n_points)
            {
              ColliderBit_error().raise(LOCAL_INFO, "ERROR! The number of grid points ("+std::to_string(n_points)+") does not agree with the number of unique 'x' and 'y' values ("+std::to_string(nx)+" and "+std::to_string(ny)+") for the interpolator '"+interp_name+"'.\n Check formatting of the file: '"+file+"'.");
            }

            // Initialize a gsl_spline pointer and two gsl_interp_accel pointers and store them
            gsl_interp_accel* x_acc = gsl_interp_accel_alloc();
            x_accels.push_back(x_acc);
            gsl_interp_accel* y_acc = gsl_interp_accel_alloc();
            y_accels.push_back(y_acc);

            gsl_spline2d* sp = gsl_spline2d_alloc(gsl_interp2d_bilinear, nx, ny);
            gsl_spline2d_init(sp, &x_data[0], &y_data[0], &interp_data[0], nx, ny);
            splines.push_back(sp);
          }
        }

        // Destructor
        ~interp2d_collection()
        {
          for (size_t interp_index = 0; interp_index < n_interpolators; ++interp_index)
          {
            gsl_spline2d_free(splines[interp_index]);
            gsl_interp_accel_free(x_accels[interp_index]);
            gsl_interp_accel_free(y_accels[interp_index]);
          }
        }

        // Evaluate a given interpolation
        double eval(double x, double y, size_t interp_index)
        {
          return gsl_spline2d_eval(splines[interp_index], x, y, x_accels[interp_index], y_accels[interp_index]);
        }
        
        // If there is only one interpolator, we can evaluate it without specifying the index
        double eval(double x, double y)
        {
          if (n_interpolators != 1)
          {
              ColliderBit_error().raise(LOCAL_INFO, "ERROR! This interp2d_collection instance contains more than one interpolator, so the interpolator index must be specified.");
          }
          return eval(x, y, 0);
        }

        // Check if point is inside interpolation range
        bool is_inside_range(double x, double y)
        {
          return ((x >= x_min) && (x <= x_max) && (y >= y_min) && (y <= y_max));
        }

    };


  }
}
