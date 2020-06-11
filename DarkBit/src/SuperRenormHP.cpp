///  \file
///
///  Super Renormalizable Higgs Portal DM specific module functions for DarkBit
///
///  *********************************************
///
///  Authors (add name and date if you modify):
///
///  \author Iñigo Saez Casares
///          (inigo.saez_casares@ens-paris-saclay.fr)
///  \date 2019 December
///
///  *********************************************

#include <algorithm>
#include <cmath>
#include <math.h>
#include <type_traits>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <functional>

#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_monte_plain.h>
#include <gsl/gsl_monte_vegas.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_min.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_sf_bessel.h>

#include "gambit/Elements/gambit_module_headers.hpp"
#include "gambit/Elements/spectrum_helpers.hpp"
#include "gambit/Utils/util_functions.hpp"
#include "gambit/Utils/ascii_table_reader.hpp"
#include "gambit/Utils/numerical_constants.hpp"
#include "gambit/Utils/statistics.hpp"
#include "gambit/DarkBit/DarkBit_rollcall.hpp"
#include "gambit/DarkBit/DarkBit_utils.hpp"

namespace Gambit
{
  namespace DarkBit
  {

    ////////////////////////////////////////////////////////////////////
    //                                                                //
    //      General Functions and Classes for Higgs Portal DM         //
    //                                                                //
    ////////////////////////////////////////////////////////////////////

    /*! \brief Supporting classes and functions for the Higgs Portal DM module.
    */

    //------------- Numerical constants and other useful things -------------//

    // masses
    const double Mp = Gambit::m_planck*1e9; // Planck mass [eV]
    const double me = Gambit::m_electron*1e9; // electron mass [eV]
    const double mH = 125.1*1e9; // Higgs boson mass [eV] (PDG 2019)
    const double mT = 172.9*1e9; // top quark mass [eV] (PDG 2019)

    // mathematical constants
    const double pi=Gambit::pi;

    // physical constants
    const double alphaEM = Gambit::alpha_EM; // fine structure constant
    const double alphaS = pi; // strong coupling constant, don't use this (take from spectrum)
    const double v = 246e9; // electroweak vev [eV]
    const double hbar_GeV = Gambit::hbar; // reduced Planck constant [GeV.s]
    const double hbar_eV = Gambit::hbar*1e9; // reduced Planck constant [eV.s]
    const double hbar_cgs = 1.054571818e-27; // reduced Planck constant [cgs]
    const double cs = Gambit::s2cm; // speed of light [cm/s]
    const double kb = Gambit::K2eV; // Boltzmann constant [eV/K]
    const double G_cgs = 6.674e-8; // Gravitational constant [cm³/g/s²]
    const double G_SI = 6.674e-11; // Gravitational constant [m³/kg/s²]

    // cosmological constants
    const double s0(2891); // current entropy density [1/cm³]
    const double rhoC(4.84e3); // current critical density [eV/cm³]

    // astrophysical constants
    const double r0(26.2225e21); // Sun's distance from the galactic center [cm]
    const double Rsun(5.9598e10); // Solar radius [cm]
    const double L0 = 2.388672e45; // solar photon luminosity [eV/s]
    const double rho0(0.3e9); // local DM density (Milky Way) [eV/cm³]

    // other constants
    const double C(50./27.); // loop function from the decay of the Higgs boson into two photons

    // Minimum finite result returnable from log(double x);
    const double logmin = log(std::numeric_limits<double>::min());


    ////////////////////////////////////////////////////////////////////
    //         Support class to handle X-ray experiments              //
    ////////////////////////////////////////////////////////////////////

    //------------- Class declaration -------------//

    class Xray
    {
      public:

        Xray(std::string experiment, double J_factor);
        double solidAngle(std::vector<double> lRange, std::vector<double> bRange);
        void set_deltaOmega();
        double getDeltaOmega() const;
        double getJ() const;
        double getEmin() const;
        double getEmax() const;
        double getDeltaE() const;
        int getFluxOrigin() const;

        double flux(double const& E);
        double sigma(double const& E);
        double fluxIntegrated(double const& E);
        double sigmaIntegrated(double const& E);
        double deltaE(double const& E);
        ~Xray();

      protected:

        double m_J; // astrophysical factor for the predicted photon flux from decaying DM
        double m_Emin; // minimum energy of the observations
        double m_Emax; // maximum energy of the observations
        double m_deltaOmega; // total solid angle of observation
        double m_deltaE; // energy resolution in percentage of the energy scale
        int m_fluxOrigin; // origin of observed flux: galactic (1), extra-galactic (2) or both (3)
        std::vector<std::vector<double> > m_lRange; // observation region in galactic coordinates (degrees)
        std::vector<std::vector<double> > m_bRange;
        std::string m_experiment;
        std::map<std::string, int> m_experimentMap;
    };

    // constructor
    Xray::Xray(std::string experiment, double J_factor) : m_experiment(experiment), m_experimentMap({{"INTEGRAL", 1}, {"HEAO", 2}}), m_J(J_factor)
    {
      switch(m_experimentMap[m_experiment])
      {
        case 1 :
          m_Emin = 20e3;
          m_Emax = 2e6;
          m_lRange = { {-30., 30.} };
          m_bRange = { {-15., 15.} };
          m_deltaE = 8e3;
          m_deltaOmega = 0.542068;
          m_fluxOrigin = 1;
          break;

        case 2 :
          m_Emin = 3e3;
          m_Emax = 60e3;
          m_lRange = { {58., 109.}, {238., 289.} };
          m_bRange = { {-90., -20.}, {20., 90.} };
          m_deltaE = 0.3;
          m_deltaOmega = 1.17135;
          m_fluxOrigin = 3;
          break;

        default :
          throw std::runtime_error("Wrong experiment name in Xray object");
          break;
      }
      //set_deltaOmega();
    }

    //------------- Function returning the energy dispersion of the instrument -------------//

    double Xray::deltaE (double const& E)
    {
      switch(m_experimentMap[m_experiment])
      {
        case 1 :
          return m_deltaE;
          break;

        case 2 :
          return m_deltaE*E;
          break;

        default :
          return 1.;
          break;
      }
    }

    // ----------- Functions to compute the solid angle of observation -------------//

    // auxiliary function for gsl integration
    double deltaOmega (double x[], size_t dim, void *p)
    {
      (void)(p);
      (void)(dim);
      return cos(x[1]);
    }

    // computes the solide angle for a given galactic coordinates range (in degrees)
    double Xray::solidAngle(std::vector<double> lRange, std::vector<double> bRange)
    {
      const size_t dim = 2, calls = 1e8;
      const double xl[dim] = {lRange[0]*pi/180., bRange[0]*pi/180.}, xu[dim] = {lRange[1]*pi/180., bRange[1]*pi/180.};
      double result, abserr;

      gsl_monte_plain_state *s = gsl_monte_plain_alloc(dim);
      gsl_monte_plain_init(s);
      gsl_rng *r = gsl_rng_alloc(gsl_rng_taus2);

      gsl_monte_function F;

      F.f = &deltaOmega;
      F.dim = dim;
      F.params = 0;

      gsl_monte_plain_integrate(&F, xl, xu, dim, calls, r, s, &result, &abserr);

      gsl_monte_plain_free(s);  

      return result;
    }

    // sets the total solid angle of observation for the experiment
    void Xray::set_deltaOmega()
    {
      double result(0);
      for (size_t i=0; i<m_lRange.size(); ++i)
      {
        result += solidAngle(m_lRange[i], m_bRange[i]);
      }
      m_deltaOmega = result;
    }

    //------------- Functions to compute the photon flux and its standard deviation -------------//

    // differential photon flux [photons/keV/cm²/s]
    double Xray::flux(double const& E)
    {
      switch(m_experimentMap[m_experiment])
      {
        case 1 :
          /* return 4.8e-8*pow(E/100e3,-1.55) + 6.6e-8*exp(-(E-50e3)/7.5e3); */
          return 1.6e-7*exp(-(E-50e3)/7.7e3) + 0.92e-7*pow(E/100e3, -1.79) + 0.34e-7*pow(E/100e3, -0.95)*exp(-(E-100e3)/3411e3) + 67.3e-7/E;
          break;

        case 2 :
          return 7.877*pow(10., 0.87)*pow(E, -0.29)*exp(-(E/41.13e3))/E*m_deltaOmega;
          break;

        default :
          return 0.;
          break;
      }
    }

    // auxiliary function for gsl integration
    double flux_gsl(double x, void *p)
    {
      Xray *experiment = static_cast<Xray*>(p);
      return experiment->flux(x);
    }

    const double int_factor(1.1); // integration range = int_factor*energy dispersion instrument

    double Xray::fluxIntegrated(double const& E)
    {
      size_t n = 1e4;

      gsl_integration_workspace *w = gsl_integration_workspace_alloc(n);

      double epsabs = 0.;
      double epsrel = 1e-2;
      size_t limit = 1e3;
      double result, abserr;
      int key = 6;
      double delta = int_factor*deltaE(E);

      gsl_function F;
      F.function = &flux_gsl;
      F.params = this;

      gsl_integration_qag(&F, E-delta/2., E+delta/2., epsabs, epsrel, limit, key, w, &result, &abserr);

      gsl_integration_workspace_free(w);

      return result;
    }

    // standard deviation of the differential photon flux [photons/keV/cm²/s]
    double Xray::sigma(double const& E)
    {

      switch(m_experimentMap[m_experiment])
      {
        case 1 :
          /* return sqrt(pow(pow(E/100e3,-1.55),2)*pow(0.6e-5,2) + pow(4.8e-8*1.55*pow(E/100e3,-2.55),2)*pow(0.25,2) + exp(-2*(E-50e3)/7.5e3)*pow(0.5e-8, 2.) + pow(6.6e-8, 2.)*pow((E-50e3)/pow(7.5e3, 2.), 2.)*exp(-2*(E-50e3)/7.5e3)*pow(1e3, 2.)); */
          return sqrt( pow(14.6e-4/E, 2) + pow(1.6e-7*(E-50e3)/7.7e3, 2) * pow(0.7e3, 2) * exp(-2.*(E-50e3)/7.7e3) + pow(0.4e-7, 2) * exp(-2.*(E-50e3)/7.7e3) + pow(0.34e-7*pow(E/100e3, -0.95), 2) * pow((E-100e3)/3411e3, 2) * exp(-(E-100e3)/3411e3) * pow(2371e3, 2));
          break;

        case 2 :
          return (flux(E)*E/m_deltaOmega)*0.02*m_deltaOmega/E;
          break;

        default :
          return 0.;
          break;
      }
    }

    // auxiliary function for gsl integration
    double sigma_gsl (double x, void *p)
    {
      Xray *experiment = static_cast<Xray*>(p);
      return pow(experiment->sigma(x), 2.);
    }


    double Xray::sigmaIntegrated(double const& E)
    {
      size_t n = 1e4;

      gsl_integration_workspace *w = gsl_integration_workspace_alloc(n);

      double epsabs = 0.;
      double epsrel = 1e-2;
      size_t limit = 1e3;
      double result, abserr;
      int key = 6;
      double delta = int_factor*deltaE(E);

      gsl_function F;
      F.function = &sigma_gsl;
      F.params = this;

      gsl_integration_qag(&F, E-delta/2., E+delta/2., epsabs, epsrel, limit, key, w, &result, &abserr);

      gsl_integration_workspace_free(w);

      return result;
    }

    //------------- Elevator functions -------------// 

    double Xray::getDeltaOmega() const { return m_deltaOmega; }

    double Xray::getJ() const { return m_J; }

    double Xray::getEmin() const { return m_Emin; }

    double Xray::getEmax() const { return m_Emax; }

    double Xray::getDeltaE() const { return m_deltaE; }

    int Xray::getFluxOrigin() const { return m_fluxOrigin; }

    // destructor
    Xray::~Xray() { }


    ////////////////////////////////////////////////////////////////////
    //            Support class to handle Solar models                //
    ////////////////////////////////////////////////////////////////////

    // gsl error handler
    void handler_Ls (const char * reason, const char * file, int line, int gsl_errno)
    {
      if (gsl_errno == 15) 
      {
        throw gsl_errno;
      }

      if (gsl_errno == 1)
      {
        std::cerr << "gsl: " << file << ":" << line << ": ERROR: " << reason << " and  gsl_errno = " << gsl_errno << std::endl;
        std::cerr << "it seems like you are trying to scan DM masses <= 0.5eV, solar neutrino likelihoods only work for DM masses > mMin = 0.5eV!" << std::endl;
        std::cerr << "if you want to change the value of mMin, you can do it from DarkBit/src/SuperRenormHP.cpp and recompute the interpolation tables(delete the existing one and it will be done automatically), but be careful to check to convergence with the new values! you may have to increase the number of points in the interpolation routine" << std::endl;
      }
      else { std::cerr << "gsl: " << file << ":" << line << ": ERROR: " << reason << " and  gsl_errno = " << gsl_errno << std::endl; abort(); }
    }

    const double da = 1.66053906660e-24; // dalton to g 


    class StellarModel
    {
      public :

        StellarModel (std::string datafile);

        double sigmaL (double const& w, double const& r);
        double getQuantity (std::string const& quantity, double const& r);
        double L_integrated (double const& mS);
        void Ls_interpolate ();
        double Ls (double const& mS, double const& theta);
        double PhiB8 (double const& mS, double const& theta);

        ~StellarModel();

      private :

        int m_nbins;
        const std::vector<std::string> m_names = {"Mass", "Radius", "Temp", "Rho", "Pres", "Lumi", "H1", "He4", "He3", "C12", "C13", "N14", "N15", "O16", "O17", "O18", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni"};
        const std::vector<std::string> m_elements = {"H1", "He4", "He3", "C12", "C13", "N14", "N15", "O16", "O17", "O18", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni"};
        // Ionisation of species i assuming full ionisation.
        const std::vector<double> m_Z = {1.0, 2.0, 2.0, 6.0, 6.0, 7.0, 7.0, 8.0, 8.0, 8.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0, 21.0, 22.0, 23.0, 24.0, 25.0, 26.0, 27.0, 28.0};
        // Atomic weight of species i (exact weight if isotope is known OR estimate from average solar abundance from data if available OR estimate from natural terrestrial abundance).
        const std::vector<double> m_A = {1.007825, 4.002603, 3.016029, 12.000000, 13.003355, 14.003074, 15.000109, 15.994915, 16.999132, 17.999160, 20.1312812, 22.989769, 24.3055, 26.9815385, 28.085, 30.973762, 32.0675, 35.4515, 36.275403, 39.0983, 40.078, 44.955908, 47.867, 50.9415, 51.9961, 54.938044, 55.845, 58.933194, 58.6934};

        ASCIItableReader m_model;
        std::map<std::string, gsl_spline*> m_interp;
        std::map<std::string, gsl_interp_accel*> m_acc;

        gsl_spline * m_Ls_interp;
        gsl_interp_accel * m_Ls_accel;

        const std::vector<std::string> m_quantities = {"Temp", "wp", "ne", "SumNz"};

    };

    StellarModel::StellarModel (std::string datafile)
    {
      m_model = ASCIItableReader(datafile);
      m_model.setcolnames(m_names);

      m_nbins = m_model["Mass"].size();

      std::vector<double> T, wp, ne, SumNz;

      double sumz, temp, sumz2;

      for (int i=0; i<m_nbins; ++i)
      { 

        sumz = 0;
        sumz2 = 0;

        for (size_t j=0; j<m_elements.size(); ++j)
        {
          sumz += m_Z[j]*m_model[m_elements[j]][i]/m_A[j]/da*m_model["Rho"][i];
          sumz2 += m_Z[j]*m_Z[j]*m_model[m_elements[j]][i]/m_A[j]/da*m_model["Rho"][i];
        }

        temp = m_model["Temp"][i]*kb;

        T.push_back(temp);
        wp.push_back(sqrt(4*pi*alphaEM*sumz*cs*cs*cs*hbar_eV*hbar_eV*hbar_eV/me));
        ne.push_back(sumz);
        SumNz.push_back(sumz2);

      }

      std::map<std::string, std::vector<double>> quantities;

      quantities["Temp"] = T;
      quantities["wp"] = wp;
      quantities["ne"] = ne;
      quantities["SumNz"] = SumNz;

      const double * rr = m_model["Radius"].data();

      for (auto it=m_quantities.begin(); it!=m_quantities.end(); ++it)
      {
        m_interp[*it] = gsl_spline_alloc(gsl_interp_cspline, m_nbins);
        m_acc[*it] = gsl_interp_accel_alloc();
        gsl_spline_init(m_interp[*it], rr, quantities[*it].data(), m_nbins);
      }

      Ls_interpolate();
    }

    // function returning the interpolated quantities (Temp, wp, ne SumNz) at a given radius r inside the Sun
    double StellarModel::getQuantity (std::string const& quantity, double const& r)
    {
      return gsl_spline_eval(m_interp[quantity], r, m_acc[quantity]);
    }

    // function to compute the damping rate of longitudinal photons inside the Sun (inverse bremsstrahlung)
    double StellarModel::sigmaL (double const& w, double const& r)
    {
      double result = 0;

      double temp = getQuantity("Temp", r), ne = getQuantity("ne", r), wp = getQuantity("wp", r);
      double D = 64*pow(pi, 2)*pow(alphaEM, 3)*ne*getQuantity("SumNz", r);
      double N = 3*sqrt(2*pi*temp)*pow(me, 3./2.)*pow(w, 3);

      double x = w/2./temp;

      // tries to use the gsl_sf_bessel_K0 and sinh functions (fails if x is too big)
      try
      {
        double F = gsl_sf_bessel_K0(x)*sinh(x);
        result += D/N*F*pow(cs, 6)*pow(hbar_eV, 6);
      }

      // uses a first order development at high x instead of the full functions
      catch (int gsl_errno)
      {
        double F = sqrt(pi/2./x)*(1-(1/8./x));
        result += D/N*F*pow(cs, 6)*pow(hbar_eV, 6);
      }

      if (w > wp) { result += 8*pi*pow(alphaEM, 2)*ne*sqrt(1-pow(wp/w, 2))/pow(me, 3)/3.*pow(cs, 3)*pow(hbar_eV, 3); }

      return result;
    }

    StellarModel::~StellarModel() 
    {
      for (auto it=m_quantities.begin(); it!=m_quantities.end(); ++it)
      {
        gsl_spline_free (m_interp[*it]);
        gsl_interp_accel_free(m_acc[*it]);
      }
    }

    struct my_f_params { double mS; StellarModel *model; double mLim; };

    // function to be integrated for mS < w < inf and 0 < r < R0
    double myF (double const& w, double const& r, my_f_params *params)
    {
      StellarModel *model = params->model;
      double mS = params->mS;
      double N1, N2, D1, D2, D3;

      N1 = pow((w*w-mS*mS), 3./2.)*w*w;
      N2 = w*w*model->sigmaL(w, r);

      D1 = pow(2*pi, 3)*alphaEM;
      D2 = exp(w/model->getQuantity("Temp", r))-1;
      D3 = pow(w*model->sigmaL(w, r), 2)+pow(w*w-pow(model->getQuantity("wp", r), 2), 2);

      return N1*N2*pow(r, 2)/D1/D2/D3;
    }

    // gsl integrand for mS < w < mLim
    double integrand1 (double x[], size_t dim, void *p)
    {
      struct my_f_params *params = (struct my_f_params *)p;
      (void)(dim);

      double w = x[0], r = x[1];

      return myF(w, r, params);
    }

    // gsl integrand for mLim < w < inf, rescaled by a change of variables to 0 < t < 1
    double integrand2 (double x[], size_t dim, void *p)
    {
      struct my_f_params *params = (struct my_f_params *)p;
      (void)(dim);

      double mLim = params->mLim;

      double t = x[0], r = x[1];

      return myF(mLim + (1-t)/t, r, params)/pow(t, 2);
    }

    const double mSmax = 1e5; // maximum mass up to which Ls is computed, for higher masses Ls is set to zero manually in the capability function

    double StellarModel::L_integrated (double const& mS)
    {
      const size_t dim = 2, calls = 1e6;
      const double mLim = sqrt(mS*mSmax*1e1);
      double xl1[dim] = {mS, 0.0006}, xu1[dim] = {mLim, 0.9995}, xl2[dim] = {0, 0.0006}, xu2[dim] = {1, 0.9995};
      double result1, abserr1, result2, abserr2;

      gsl_monte_vegas_state *s1 = gsl_monte_vegas_alloc(dim), *s2 = gsl_monte_vegas_alloc(dim);
      gsl_monte_vegas_init(s1);
      gsl_monte_vegas_init(s2);

      gsl_rng *r = gsl_rng_alloc(gsl_rng_taus2);

      gsl_monte_function F1, F2;

      my_f_params params = {mS, this, mLim};

      F1.f = &integrand1;
      F1.dim = dim;
      F1.params = &params;

      F2.f = &integrand2;
      F2.dim = dim;
      F2.params = &params;

      gsl_monte_vegas_integrate(&F1, xl1, xu1, dim, calls, r, s1, &result1, &abserr1);
      gsl_monte_vegas_integrate(&F2, xl2, xu2, dim, calls, r, s2, &result2, &abserr2);

      std::cout << mS << " " << result1 << " " << abserr1/result1 << " " << result2 << " " << abserr2/result2 << std::endl;

      gsl_monte_vegas_free(s1);
      gsl_monte_vegas_free(s2);

      return 4*pi*pow(Rsun, 3)*(result1+result2)/pow(cs, 3)/pow(hbar_eV, 4);
    }

    void StellarModel::Ls_interpolate ()
    {
      const int nPoints = 100;
      const double mMin = 1e-9, mMax = mSmax;
      const double deltaM = log10(mMax/mMin)/nPoints;
      const std::string filename = GAMBIT_DIR "/DarkBit/data/SuperRenormHP_Ls.dat";
      std::vector<double> mS, Ls;

      // checks if an interpolation table file for Ls already exists and reads it if that is the case
      if (Utils::file_exists(filename))
      {
        ASCIItableReader table = ASCIItableReader(filename);
        table.setcolnames({"mS", "Ls"});
        mS = table["mS"];
        Ls = table["Ls"];
      }

      // builds an interpolation table for Ls and writes it in a file
      else
      {
        for (int i=0; i<=nPoints; ++i)
        {
          mS.push_back(mMin*pow(10, deltaM*i));
        }

        for (auto it=mS.begin(); it!=mS.end(); ++it)
        {
          Ls.push_back(L_integrated(*it));
        }

        std::ofstream fout (filename.c_str());

        for (size_t i=0; i<mS.size(); ++i)
        {
          fout << mS[i] << "  " << Ls[i] << std::endl;
        }

        fout.clear(); fout.close();
      }

      // interpolates Ls from the interpolation table
      m_Ls_interp = gsl_spline_alloc(gsl_interp_cspline, nPoints+1);
      m_Ls_accel = gsl_interp_accel_alloc();
      gsl_spline_init(m_Ls_interp, mS.data(), Ls.data(), nPoints+1);
    }

    // returns the value of the interpolated Ls for a given set of parameters (mS, theta)
    double StellarModel::Ls (double const& mS, double const& theta)
    {
      return gsl_spline_eval(m_Ls_interp, mS, m_Ls_accel)*pow(me/v*theta, 2);
    }


    ////////////////////////////////////////////////////////////////////
    //                   Capability functions                         //
    ////////////////////////////////////////////////////////////////////

    //------------- Process catalogue -------------//

    void TH_ProcessCatalog_SuperRenormHP(DarkBit::TH_ProcessCatalog &result)
    {
      using namespace Pipes::TH_ProcessCatalog_SuperRenormHP;
      using std::vector;
      using std::string;

      // Initialize empty catalog and decay channel
      TH_ProcessCatalog catalog;
      TH_Process process_dec("S");

      ///////////////////////////////////////
      // Import particle masses and couplings
      ///////////////////////////////////////

      // Convenience macros
      #define getSMmass(Name, spinX2)                                           \
       catalog.particleProperties.insert(std::pair<string, TH_ParticleProperty> \
       (Name , TH_ParticleProperty(SM.get(Par::Pole_Mass,Name), spinX2)));
      #define addParticle(Name, Mass, spinX2)                                   \
       catalog.particleProperties.insert(std::pair<string, TH_ParticleProperty> \
       (Name , TH_ParticleProperty(Mass, spinX2)));

      // Import Spectrum objects
      const Spectrum& spec  = *Dep::SuperRenormHP_spectrum;
      const SubSpectrum& he = spec.get_HE();
      const SubSpectrum& SM = spec.get_LE();
      const SMInputs& SMI   = spec.get_SMInputs();

      // Import couplings
      double theta = he.get(Par::dimensionless, "theta");
      double vev = he.get(Par::mass1,"vev");
      double C = 50./27.; // loop factor for the decay into two photons
      double alpha = 1./SMI.alphainv; // alpha_EM(mZ)^MSbar (5 active flavours)

      // Get SM pole masses
      getSMmass("e-_1",     1)
      getSMmass("e+_1",     1)
      getSMmass("e-_2",     1)
      getSMmass("e+_2",     1)
      getSMmass("e-_3",     1)
      getSMmass("e+_3",     1)
      getSMmass("Z0",     2)
      getSMmass("W+",     2)
      getSMmass("W-",     2)
      getSMmass("g",      2)
      getSMmass("gamma",  2)
      getSMmass("u_3",      1)
      getSMmass("ubar_3",   1)
      getSMmass("d_3",      1)
      getSMmass("dbar_3",   1)

      // Pole masses not available for the light quarks.
      addParticle("u_1"   , SMI.mU,  1) // mu(2 GeV)^MS-bar, not pole mass
      addParticle("ubar_1", SMI.mU,  1) // mu(2 GeV)^MS-bar, not pole mass
      addParticle("d_1"   , SMI.mD,  1) // md(2 GeV)^MS-bar, not pole mass
      addParticle("dbar_1", SMI.mD,  1) // md(2 GeV)^MS-bar, not pole mass
      addParticle("u_2"   , SMI.mCmC,1) // mc(mc)^MS-bar, not pole mass
      addParticle("ubar_2", SMI.mCmC,1) // mc(mc)^MS-bar, not pole mass
      addParticle("d_2"   , SMI.mS,  1) // ms(2 GeV)^MS-bar, not pole mass
      addParticle("dbar_2", SMI.mS,  1) // ms(2 GeV)^MS-bar, not pole mass
      /* double alpha_s = SMI.alphaS;      // alpha_s(mZ)^MSbar */

      // Masses for neutrino flavour eigenstates. Set to zero.
      // (presently not required)
      addParticle("nu_e",     0.0, 1)
      addParticle("nubar_e",  0.0, 1)
      addParticle("nu_mu",    0.0, 1)
      addParticle("nubar_mu", 0.0, 1)
      addParticle("nu_tau",   0.0, 1)
      addParticle("nubar_tau",0.0, 1)

      // Higgs-sector masses
      double mS = spec.get(Par::Pole_Mass,"S");
      double mH = spec.get(Par::Pole_Mass,"h0_1");
      addParticle("S",        mS, 0)  // Scalar DM
      addParticle("h0_1",     mH, 0)  // SM-like Higgs
      addParticle("pi0",   meson_masses.pi0,       0)
      addParticle("pi+",   meson_masses.pi_plus,   0)
      addParticle("pi-",   meson_masses.pi_minus,  0)
      addParticle("eta",   meson_masses.eta,       0)
      addParticle("rho0",  meson_masses.rho0,      1)
      addParticle("rho+",  meson_masses.rho_plus,  1)
      addParticle("rho-",  meson_masses.rho_minus, 1)
      addParticle("omega", meson_masses.omega,     1)

      // Get rid of convenience macros
      #undef getSMmass
      #undef addParticle


      // decay into two photons (through a loop of heavy fermions)
      double gamma = (theta*theta*alpha*alpha*mS*mS*mS*C*C)/(256.*pi*pi*pi*vev*vev);

      TH_Channel dec_channel(daFunk::vec<string>("gamma", "gamma"), daFunk::cnst(gamma));
      process_dec.channelList.push_back(dec_channel);

      //////////////////////////////
      // Import Decay information //
      //////////////////////////////

      // Import decay table from DecayBit
      const DecayTable* tbl = &(*Dep::decay_rates);

      // Set of imported decays
      std::set<string> importedDecays;

      // Minimum branching ratio
      double minBranching = 0.;

      // Import relevant decays (only Higgs and subsequent decays)
      using DarkBit_utils::ImportDecays;
      // Notes: Virtual Higgs decays into offshell W+W- final states are not
      // imported.  All other channels are correspondingly rescaled.  Decay
      // into SS final states is accounted for, leading to zero photons.
      ImportDecays("h0_1", catalog, importedDecays, tbl, minBranching, daFunk::vec<std::string>("Z0", "W+", "W-"));

      // Add the decay process to the catalog
      catalog.processList.push_back(process_dec);

      // Validate
      /* catalog.validate(); */

      result = catalog;
    } // function TH_ProcessCatalog_SuperRenormHP


    //------------- Functions to compute the age of the Universe at a given redshift -------------//

    // useful structure
    struct cosmology_params {double OmegaM; double OmegaR; double OmegaLambda; double H0;};

    // auxiliary function for gsl integration
    double age_f (double x, void *p)
    {
      cosmology_params *params = static_cast<cosmology_params*>(p);
      double OmegaM = params->OmegaM;
      double OmegaLambda = params->OmegaLambda;
      double OmegaR = params->OmegaR;
      double OmegaK = 0.;

      return 1./sqrt( OmegaM*pow(1+x, 5.) + OmegaLambda*pow(1+x, 2.) + OmegaK*pow(1+x, 4.) + OmegaR*pow(1+x, 6.) );
    }

    // computes the age of the Universe at a given redshift ([0] age, [1] abserr)
    std::vector<double> ageUniverse (double redshift, double OmegaM, double OmegaR, double OmegaLambda, double H0)
    {
      size_t n = 1e4; 

      gsl_integration_workspace *w =  gsl_integration_workspace_alloc(n);

      double epsabs = 0.;
      double epsrel = 1e-3;
      size_t limit = 1e3;
      double result, abserr;

      cosmology_params params = {OmegaM, OmegaR, OmegaLambda, H0};

      gsl_function F;
      F.function = &age_f;
      F.params = &params;

      gsl_integration_qagiu(&F, redshift, epsabs, epsrel, limit, w, &result, &abserr);

      gsl_integration_workspace_free(w);

      return {result/H0, abserr/H0};
    }

    //------------- Functions to compute X-ray likelihoods -------------//


    // useful structure
    struct XrayLikelihood_params {double mass; double gamma; double density; Xray experiment; double H0; double OmegaM; double OmegaR; double OmegaLambda; double OmegaDM;};

    // extra-galactic contribution to the differential photon flux [photons/eV/cm²/s]
    double dPhiEG(double const& E, XrayLikelihood_params *params)
    {
      double mass = params->mass, gamma = params->gamma, density = params->density;
      Xray experiment = params->experiment;

      double H0 = params->H0;
      double OmegaM = params->OmegaM;
      double OmegaR = params->OmegaR;
      double OmegaLambda = params->OmegaLambda;
      double OmegaK = 0.;

      double x = mass/2./E;
      double z = x - 1.;

      double t = ageUniverse(z, OmegaM, OmegaR, OmegaLambda, H0)[0];

      return experiment.getDeltaOmega()*2.*1./(4*pi)*(gamma*density*cs*exp(-gamma*t))/(mass*H0*E)/sqrt( OmegaM*pow(x, 3.) + OmegaLambda + OmegaK*pow(x, 2.) + OmegaR*pow(x, 4.) );
    }

    const double s(1./3.); // standard deviation of the gaussian for the galactic emission line = s*energy dispersion instrument
    // galactic (Milky Way) contribution to the differential photon flux [photons/eV/cm²/s]
    double dPhiG(double const& E, XrayLikelihood_params *params)
    {
      double mass = params->mass, gamma = params->gamma, density = params->density;
      Xray experiment = params->experiment;

      double J = experiment.getJ();

      double H0 = params->H0;
      double OmegaM = params->OmegaM;
      double OmegaR = params->OmegaR;
      double OmegaLambda = params->OmegaLambda;

      double t0 = ageUniverse(0., OmegaM, OmegaR, OmegaLambda, H0)[0];
      double OmegaDM = params->OmegaDM;

      double sigma = s*experiment.deltaE(E); // standard deviation of the gaussian modelling the enery dispersion of the instrument

      return 2.*(gamma*J*density*exp(-t0*gamma))/(4.*pi*mass*OmegaDM*rhoC)/sqrt(2*pi*sigma*sigma)*exp(-pow(E-mass/2.,2)/(2*sigma*sigma));
    }

    // total predicted differential photon flux for a given X-ray experiment [photons/eV/cm²/s]
    double XrayPrediction(double const& E, XrayLikelihood_params *params)
    {
      Xray experiment = params->experiment;
      switch(experiment.getFluxOrigin())
      {
        case 1 :
          return dPhiG(E, params);

        case 2 :
          return dPhiEG(E, params);

        case 3 :
          return dPhiG(E, params) + dPhiEG(E, params);

        default :
          throw std::runtime_error("Wrong value for m_fluxOrigin in Xray class, allowed values are 1 (galactic flux), 2 (extra-galactic flux) and 3 (both)");
          break;
      }

    }

    // auxiliary function for gsl integration
    double XrayPredictionIntegrated_gsl(double x, void *p)
    {
      XrayLikelihood_params *params = static_cast<XrayLikelihood_params*>(p);

      return XrayPrediction(x, params);
    }

    // total predicted photon flux integrated over an interval deltaE, centered around E [photons/cm²/s]
    double XrayPredictionIntegrated(double const& E, XrayLikelihood_params *params)
    {
      size_t n = 1e4;

      gsl_integration_workspace *w =  gsl_integration_workspace_alloc(n);

      double epsabs = 0.;
      double epsrel = 1e-2;
      size_t limit = 1e3;
      double result, abserr;
      int key = 6;
      Xray experiment = params->experiment;
      double delta = int_factor*experiment.deltaE(E);

      gsl_function F;
      F.function = &XrayPredictionIntegrated_gsl;
      F.params = params;

      gsl_integration_qag(&F, E-delta/2., E+delta/2., epsabs, epsrel, limit, key, w, &result, &abserr);

      gsl_integration_workspace_free(w);

      return result;
    }

    void SuperRenormHP_DecayFluxG (double &result)
    {
      using namespace Pipes::SuperRenormHP_DecayFluxG;

      double density = *Dep::DM_relic_density*1e9;

      double OmegaDM = *Dep::Omega0_cdm, H0 = *Dep::H0;

      double OmegaM = *Dep::Omega0_m, OmegaR = *Dep::Omega0_r;

      double OmegaLambda = *Dep::Omega0_Lambda;

      double RhoC = 3*pow(H0, 2)/(8*pi*G_SI);

      std::string DM_ID = *Dep::DarkMatter_ID;

      TH_ProcessCatalog catalog = *Dep::TH_ProcessCatalog;
      auto f = catalog.getProcess(DM_ID).find({"gamma", "gamma"})->genRate;
      auto fb = f->bind();
      double gamma = fb->eval();

      double mass = catalog.getParticleProperty(DM_ID).mass*1e9;

      double t0 = ageUniverse(0., OmegaM, OmegaR, OmegaLambda, H0)[0];

      result = 2.*(gamma*density*exp(-t0*gamma))/(4.*pi*mass*OmegaDM*RhoC);
    }

    // auxiliary function for gsl minimization returning the log-likelihood for a given X-ray experiment
    double XrayLogLikelihood(double E, void *p)
    {
      XrayLikelihood_params *params = static_cast<XrayLikelihood_params*>(p);
      Xray experiment = params->experiment;

      double data = experiment.fluxIntegrated(E);
      double sigma = experiment.sigmaIntegrated(E);
      double prediction = XrayPredictionIntegrated(E, params);

      return (prediction>=data) ? -pow(data-prediction,2.)/(2.*sigma*sigma) : 0.;
    }

    // computes the energy E which minimizes the log-likelihood
    double minimizeLogLikelihood(XrayLikelihood_params *params)
    {
      int status;
      int iter = 0, max_iter = 100;
      const gsl_min_fminimizer_type *T;
      gsl_min_fminimizer *s;
      Xray experiment = params->experiment;
      double mass = params->mass;
      double Emin = experiment.getEmin(), Emax = experiment.getEmax();
      double a = Emin+experiment.deltaE(Emin), b = fmin(mass/2., Emax-experiment.deltaE(Emax));
      double m = (a+b)/2.;
      gsl_function F;

      F.function = &XrayLogLikelihood;
      F.params = params;
      T = gsl_min_fminimizer_brent;
      s = gsl_min_fminimizer_alloc(T);
      gsl_min_fminimizer_set (s, &F, m, a, b);

      do
      {
        iter++;
        status = gsl_min_fminimizer_iterate (s);

        m = gsl_min_fminimizer_x_minimum (s);
        a = gsl_min_fminimizer_x_lower (s);
        b = gsl_min_fminimizer_x_upper (s);

        status = gsl_min_test_interval (a, b, 0.1, 0);
      } while (status == GSL_CONTINUE && iter < max_iter);

      gsl_min_fminimizer_free (s);

      return m;
    }

    // capability function to provide the initial energy density as produced by the freeze-in mechanism
    void SuperRenormHP_relic_density (double &result)
    {
      using namespace Pipes::SuperRenormHP_relic_density;
      double mS = *Param["mS"], theta = *Param["theta"];
      result = 4e11*theta*theta*s0*mS;
    }

    // Linear interpolation in lin-log space.
    double interpolate(double x, const std::vector<double> & xlist,
            const std::vector<double> & ylist, bool zerobound)
    {
        double x0, x1, y0, y1;
        int i = 1;
        if (zerobound)
        {
            if (x<xlist.front()) return 0;
            if (x>xlist.back()) return 0;
        }
        else
        {
            if (x<xlist.front()) return ylist.front();
            if (x>xlist.back()) return ylist.back();
        }
        // Find min i such that xlist[i]>=x.
        for (; xlist[i] < x; i++) {};
        x0 = xlist[i-1];
        x1 = xlist[i];
        y0 = ylist[i-1];
        y1 = ylist[i];
        // lin-vs-log interpolation for lnL vs flux
        return y0 + (y1-y0) * log(x/x0) / log(x1/x0);
    }

    void get_J_factor_INTEGRAL_CO (double &result)
    {
      using namespace Pipes::get_J_factor_INTEGRAL_CO;

      GalacticHaloProperties halo = *Dep::GalacticHalo;

      daFunk::Funk profile = halo.DensityProfile;

      std::vector<double> rho;
      auto r = daFunk::logspace(-3, 2, 100);
      double r_sun = halo.r_sun;

      for ( size_t i = 0; i<r.size(); i++ )
      {
        rho.push_back(profile->bind("r")->eval(r[i]));
      }

      std::vector<double> phi_pre;
      std::vector<double> intensity;

      BEreq::los_integral(byVal(r), byVal(rho), byVal(r_sun), phi_pre, intensity);

      auto emission = std::pair< std::vector<double>, std::vector<double> > (phi_pre, intensity);

      ASCIItableReader ROI = ASCIItableReader(GAMBIT_DIR "/DarkBit/data/INTEGRAL/ROI_CO.txt");
      ROI.setcolnames({"phi", "weight"});
      std::vector<double> phi = ROI["phi"], weight = ROI["weight"];

      double J = 0;
      for ( size_t i = 0; i < phi.size(); i++ )
      {
        J += interpolate(phi[i], emission.first, emission.second, true)*weight[i];
      }

      result = J;
    }

    // gsl error handler
    void handler (const char * reason, const char * file, int line, int gsl_errno)
    {
      if (gsl_errno == 4)
      {
        throw gsl_errno;
      }
      else { std::cerr << "gsl: " << file << ":" << line << ": ERROR: " << reason << std::endl; abort(); }
    }

    // capability function to compute the X-ray log-likelihood from the INTEGRAL experiment
    void calc_lnL_INTEGRAL_CO(double &result)
    {
      using namespace Pipes::calc_lnL_INTEGRAL_CO;

      double J_factor = *Dep::J_factor_INTEGRAL_CO;

      static Xray experiment = Xray("INTEGRAL", J_factor);

      double density = *Dep::DM_relic_density*1e9;

      double OmegaDM = *Dep::Omega0_cdm, H0 = *Dep::H0;

      double OmegaM = *Dep::Omega0_m, OmegaR = *Dep::Omega0_r;

      double OmegaLambda = *Dep::Omega0_Lambda;

      std::string DM_ID = *Dep::DarkMatter_ID;

      TH_ProcessCatalog catalog = *Dep::TH_ProcessCatalog;
      auto f = catalog.getProcess(DM_ID).find({"gamma", "gamma"})->genRate;
      auto fb = f->bind();
      double gamma = fb->eval();

      double mass = catalog.getParticleProperty(DM_ID).mass*1e9;

      XrayLikelihood_params params = {mass, gamma/hbar_GeV, density, experiment, H0*1e-19/3.085, OmegaM, OmegaR, OmegaLambda, OmegaDM};

      double Emin = experiment.getEmin(), Emax = experiment.getEmax(), E, lik1, lik2;

      // no constraints available above the electron threshold, we need to take into account the decay into charged particles and their subsequent FSR
      if (mass >= 1e6) { result = 0; }

      else if (mass >= 2.*Emin)
      {
        // modifies the gsl error handler and stores the default one
        gsl_error_handler_t *old_handler = gsl_set_error_handler (&handler);
        try
        {
          E = minimizeLogLikelihood(&params);
          result = XrayLogLikelihood(E, &params);
        }

        catch (int gsl_errno)
        {
          lik1 = XrayLogLikelihood(Emin+experiment.deltaE(Emin), &params);
          lik2 = XrayLogLikelihood(fmin(mass/2., Emax-experiment.deltaE(Emax)), &params);
          result = fmin(lik1, lik2);
        }
        // restores the default gsl error handler
        gsl_set_error_handler (old_handler);
      }

      else { result = 0; }
    }

    void get_J_factor_INTEGRAL_ang_b (std::vector<double> &result)
    {
      using namespace Pipes::get_J_factor_INTEGRAL_ang_b;

      GalacticHaloProperties halo = *Dep::GalacticHalo;

      daFunk::Funk profile = halo.DensityProfile;

      std::vector<double> rho;
      auto r = daFunk::logspace(-3, 2, 100);
      double r_sun = halo.r_sun;

      for ( size_t i = 0; i<r.size(); i++ )
      {
        rho.push_back(profile->bind("r")->eval(r[i]));
      }

      std::vector<double> phi_pre;
      std::vector<double> intensity;

      BEreq::los_integral(byVal(r), byVal(rho), byVal(r_sun), phi_pre, intensity);

      auto emission = std::pair< std::vector<double>, std::vector<double> > (phi_pre, intensity);

      ASCIItableReader ROI_1 = ASCIItableReader(GAMBIT_DIR "/DarkBit/data/INTEGRAL/ROI_ang_b_1.txt");
      ASCIItableReader ROI_2 = ASCIItableReader(GAMBIT_DIR "/DarkBit/data/INTEGRAL/ROI_ang_b_2.txt");

      ROI_1.setcolnames({"phi", "weight"});
      ROI_2.setcolnames({"phi", "weight"});

      std::vector<double> phi_1 = ROI_1["phi"], weight_1 = ROI_1["weight"];
      std::vector<double> phi_2 = ROI_2["phi"], weight_2 = ROI_2["weight"];

      double J_1 = 0, J_2 = 0;

      for ( size_t i = 0; i < phi_1.size(); i++ )
      {
        J_1 += interpolate(phi_1[i], emission.first, emission.second, true)*weight_1[i];
      }

      for ( size_t i = 0; i < phi_2.size(); i++ )
      {
        J_2 += interpolate(phi_2[i], emission.first, emission.second, true)*weight_2[i];
      }

      result = {J_1, J_2};
    }

    void calc_lnL_INTEGRAL_ang_b (double &result)
    {
      using namespace Pipes::calc_lnL_INTEGRAL_ang_b;

      std::vector<double> J_factor = *Dep::J_factor_INTEGRAL_ang_b;

      double mass = *Dep::DM_mass;

      double DecayFluxG = *Dep::DM_DecayFluxG;

      static ASCIItableReader INTEGRAL = ASCIItableReader(GAMBIT_DIR "/DarkBit/data/INTEGRAL/INTEGRAL_b.dat");

      INTEGRAL.setcolnames({"Emin", "Emax", "Flux", "Sigma"});

      static std::vector<double> Emin = INTEGRAL["Emin"], Emax = INTEGRAL["Emax"], Flux = INTEGRAL["Flux"], Sigma = INTEGRAL["Sigma"];

      std::vector<double> Omega = {1.6119, 4.1858};

      // no constraints available above the electron threshold, we need to take into account the decay into charged particles and their subsequent FSR
      if (mass >= 1e6) { result = 0; }

      else if (mass < 2.**std::min_element(Emin.begin(), Emin.end())) { result = 0; }

      else
      {
        std::vector<double> likelihood;

        double PredictedFlux, ObservedFlux, Error;

        for (size_t i = 0; i < Emin.size()-1; ++i)
        {
          PredictedFlux = ( (mass >= 2*Emin[i]) && (mass < 2*Emax[i]) ) ? DecayFluxG*J_factor[0]/Omega[0] : 0;
          ObservedFlux = Flux[i];
          Error = Sigma[i];
          likelihood.push_back( (PredictedFlux < ObservedFlux) ? 1 : exp(-pow(ObservedFlux - PredictedFlux, 2)/pow(Error, 2)) );
        }

        PredictedFlux = ( (mass >= 2*Emin.back()) && (mass < 2*Emax.back()) ) ? DecayFluxG*J_factor[1]/Omega[1] : 0;
        ObservedFlux = Flux.back();
        Error = Sigma.back();
        likelihood.push_back( (PredictedFlux < ObservedFlux) ? 1 : exp(-pow(ObservedFlux - PredictedFlux, 2)/pow(Error, 2)) );

        result = log(*std::min_element(likelihood.begin(), likelihood.end()));
      }

    }

    void get_J_factor_INTEGRAL_ang_l (std::vector<double> &result)
    {
      using namespace Pipes::get_J_factor_INTEGRAL_ang_l;

      GalacticHaloProperties halo = *Dep::GalacticHalo;

      daFunk::Funk profile = halo.DensityProfile;

      std::vector<double> rho;
      auto r = daFunk::logspace(-3, 2, 100);
      double r_sun = halo.r_sun;

      for ( size_t i = 0; i<r.size(); i++ )
      {
        rho.push_back(profile->bind("r")->eval(r[i]));
      }

      std::vector<double> phi_pre;
      std::vector<double> intensity;

      BEreq::los_integral(byVal(r), byVal(rho), byVal(r_sun), phi_pre, intensity);

      auto emission = std::pair< std::vector<double>, std::vector<double> > (phi_pre, intensity);

      ASCIItableReader ROI_1 = ASCIItableReader(GAMBIT_DIR "/DarkBit/data/INTEGRAL/ROI_ang_l_1.txt");
      ASCIItableReader ROI_2 = ASCIItableReader(GAMBIT_DIR "/DarkBit/data/INTEGRAL/ROI_ang_l_2.txt");

      ROI_1.setcolnames({"phi", "weight"});
      ROI_2.setcolnames({"phi", "weight"});

      std::vector<double> phi_1 = ROI_1["phi"], weight_1 = ROI_1["weight"];
      std::vector<double> phi_2 = ROI_2["phi"], weight_2 = ROI_2["weight"];

      double J_1 = 0, J_2 = 0;

      for ( size_t i = 0; i < phi_1.size(); i++ )
      {
        J_1 += interpolate(phi_1[i], emission.first, emission.second, true)*weight_1[i];
      }

      for ( size_t i = 0; i < phi_2.size(); i++ )
      {
        J_2 += interpolate(phi_2[i], emission.first, emission.second, true)*weight_2[i];
      }

      result = {J_1, J_2};
    }

    void calc_lnL_INTEGRAL_ang_l (double &result)
    {
      using namespace Pipes::calc_lnL_INTEGRAL_ang_l;

      std::vector<double> J_factor = *Dep::J_factor_INTEGRAL_ang_l;

      double mass = *Dep::DM_mass;

      double DecayFluxG = *Dep::DM_DecayFluxG;

      static ASCIItableReader INTEGRAL = ASCIItableReader(GAMBIT_DIR "/DarkBit/data/INTEGRAL/INTEGRAL_l.dat");

      INTEGRAL.setcolnames({"Emin", "Emax", "Flux", "Sigma"});

      static std::vector<double> Emin = INTEGRAL["Emin"], Emax = INTEGRAL["Emax"], Flux = INTEGRAL["Flux"], Sigma = INTEGRAL["Sigma"];

      std::vector<double> Omega = {1.4224, 1.7919};

      // no constraints available above the electron threshold, we need to take into account the decay into charged particles and their subsequent FSR
      if (mass >= 1e6) { result = 0; }

      else if (mass < 2.**std::min_element(Emin.begin(), Emin.end())) { result = 0; }

      else
      {
        std::vector<double> likelihood;

        double PredictedFlux, ObservedFlux, Error;

        for (size_t i = 0; i < Emin.size()-1; ++i)
        {
          PredictedFlux = ( (mass >= 2*Emin[i]) && (mass < 2*Emax[i]) ) ? DecayFluxG*J_factor[0]/Omega[0] : 0;
          ObservedFlux = Flux[i];
          Error = Sigma[i];
          likelihood.push_back( (PredictedFlux < ObservedFlux) ? 1 : exp(-pow(ObservedFlux - PredictedFlux, 2)/pow(Error, 2)) );
        }

        PredictedFlux = ( (mass >= 2*Emin.back()) && (mass < 2*Emax.back()) ) ? DecayFluxG*J_factor[1]/Omega[1] : 0;
        ObservedFlux = Flux.back();
        Error = Sigma.back();
        likelihood.push_back( (PredictedFlux < ObservedFlux) ? 1 : exp(-pow(ObservedFlux - PredictedFlux, 2)/pow(Error, 2)) );

        result = log(*std::min_element(likelihood.begin(), likelihood.end()));
      }

    }

    void get_J_factor_HEAO (double &result)
    {
      using namespace Pipes::get_J_factor_HEAO;

      GalacticHaloProperties halo = *Dep::GalacticHalo;

      daFunk::Funk profile = halo.DensityProfile;

      std::vector<double> rho;
      auto r = daFunk::logspace(-3, 2, 100);
      double r_sun = halo.r_sun;

      for ( size_t i = 0; i<r.size(); i++ )
      {
        rho.push_back(profile->bind("r")->eval(r[i]));
        /* rho.push_back(pow(profile->bind("r")->eval(r[i]), 2)); */
      }

      std::vector<double> phi_pre;
      std::vector<double> intensity;

      BEreq::los_integral(byVal(r), byVal(rho), byVal(r_sun), phi_pre, intensity);

      auto emission = std::pair< std::vector<double>, std::vector<double> > (phi_pre, intensity);

      ASCIItableReader ROI = ASCIItableReader(GAMBIT_DIR "/DarkBit/data/HEAO/ROI.txt");
      ROI.setcolnames({"phi", "weight"});
      std::vector<double> phi = ROI["phi"], weight = ROI["weight"];

      double J = 0;
      for ( size_t i = 0; i < phi.size(); i++ )
      {
        J += interpolate(phi[i], emission.first, emission.second, true)*weight[i];// /Omega;
      }

      double rho_sun = profile->bind("r")->eval(r_sun);

      result = J;

      /* std::cout << "J = " << result/r_sun/rho_sun << std::endl; */
    }

    // capability function to compute the X-ray log-likelihood from the HEAO-1 A2 experiment
    void calc_lnL_HEAO(double &result)
    {
      using namespace Pipes::calc_lnL_HEAO;

      static Xray experiment = Xray("HEAO", 9.894); // J in Gev kpc /cm^3

      double density = *Dep::DM_relic_density*1e9;

      double OmegaDM = *Dep::Omega0_cdm, H0 = *Dep::H0;

      double OmegaM = *Dep::Omega0_m, OmegaR = *Dep::Omega0_r;

      double OmegaLambda = *Dep::Omega0_Lambda;

      std::string DM_ID = *Dep::DarkMatter_ID;

      TH_ProcessCatalog catalog = *Dep::TH_ProcessCatalog;
      auto f = catalog.getProcess(DM_ID).find({"gamma", "gamma"})->genRate;
      auto fb = f->bind();
      double gamma = fb->eval();

      double mass = catalog.getParticleProperty(DM_ID).mass*1e9;

      XrayLikelihood_params params = {mass, gamma/hbar_GeV, density, experiment, H0*1e-19/3.085, OmegaM, OmegaR, OmegaLambda, OmegaDM};

      const double Emin = experiment.getEmin(), Emax = experiment.getEmax();
      double E, lik1, lik2;

      // no constraints available above the electron threshold, we need to take into account the decay into charged particles and their subsequent FSR
      if (mass >= 1e6) { result = 0; }

      else if (mass >= 2.*Emin)
      {
        // modifies the gsl error handler and stores the default one
        gsl_error_handler_t *old_handler = gsl_set_error_handler (&handler);
        try
        {
          E = minimizeLogLikelihood(&params);
          result = XrayLogLikelihood(E, &params);
        }

        catch (int gsl_errno)
        {
          lik1 = XrayLogLikelihood(Emin+experiment.deltaE(Emin), &params);
          lik2 = XrayLogLikelihood(fmin(mass/2., Emax-experiment.deltaE(Emax)), &params);
          result = fmin(lik1, lik2);
        }
        // restores the default gsl error handler
        gsl_set_error_handler (old_handler);
      }

      else { result = 0; }
    }

    //------------- Functions to compute stellar coolling likelihoods -------------//

    // capability function to compute the solar DM luminosity
    void SuperRenormHP_solar_luminosity (double &result)
    {
      using namespace Pipes::SuperRenormHP_solar_luminosity;
      const double mS = *Param["mS"]*1e9, theta = *Param["theta"];

      gsl_error_handler_t *old_handler = gsl_set_error_handler (&handler_Ls);

      static StellarModel model (GAMBIT_DIR "/DarkBit/data/SolarModel_AGSS09met.dat");

      if (mS < mSmax) { result = model.Ls(mS, theta); }

      else { result = 0; }

      gsl_set_error_handler (old_handler);
    }

    // capability function to compute the likelihood from the solar luminosity limit ( L_DM < 0.1*L0 ) (conservative limit)
    void calc_lnL_solar_luminosity (double &result)
    {
      using namespace Pipes::calc_lnL_solar_luminosity;

      const double limit = 0.1*L0;
      const double Ls = *Dep::solar_DM_luminosity;

      if (Ls > limit) { result = logmin; }

      else { result = 0; }
    }

    // capability function to compute the predicted solar neutrino flux (B8)
    void SuperRenormHP_solar_neutrino_flux_B8 (double &result)
    {
      using namespace Pipes::SuperRenormHP_solar_neutrino_flux_B8;

      const double Ls = *Dep::solar_DM_luminosity;
      const double alpha = runOptions->getValueOrDef<double>(4., "alpha");

      const double Phi0 = 4.95e6;

      result = Phi0*pow(1+Ls/L0, alpha);
    }

    // capability function to compute the predicted solar neutrino flux (Be7)
    void SuperRenormHP_solar_neutrino_flux_Be7 (double &result)
    {
      using namespace Pipes::SuperRenormHP_solar_neutrino_flux_Be7;

      const double Ls = *Dep::solar_DM_luminosity;
      const double alpha = runOptions->getValueOrDef<double>(4., "alpha");

      const double Phi0 = 4.71e9;

      result = Phi0*pow(1+Ls/L0, alpha);
    }

    // capability function to compute the likelihood from the solar B8 neutrino flux
    void calc_lnL_solar_neutrino_B8 (double &result)
    {
      using namespace Pipes::calc_lnL_solar_neutrino_B8;

      const bool profile = runOptions->getValueOrDef<bool>(false, "profile_systematics");
      const double Phi_predicted = *Dep::solar_neutrino_flux_B8;

      const double Phi_obs = 5e6;
      const double sigma_obs = 0.03*Phi_obs, sigma_theo = 0.14*Phi_predicted;

      result = Stats::gaussian_upper_limit(Phi_predicted, Phi_obs, sigma_theo, sigma_obs, profile);
    }

    // capability function to compute the likelihood from the solar Be7 neutrino flux
    void calc_lnL_solar_neutrino_Be7 (double &result)
    {
      using namespace Pipes::calc_lnL_solar_neutrino_Be7;

      const bool profile = runOptions->getValueOrDef<bool>(false, "profile_systematics");
      const double Phi_predicted = *Dep::solar_neutrino_flux_Be7;

      const double Phi_obs = 4.82e9;
      const double sigma_obs = 0.05*Phi_obs, sigma_theo = 0.07*Phi_predicted;

      result = Stats::gaussian_upper_limit(Phi_predicted, Phi_obs, sigma_theo, sigma_obs, profile);
    }

    //------------- Functions to compute short range forces likelihoods -------------// 

    // capability to provide the Higgs-Nucleon coupling constant fN, such as described in arXiv:1306.4710
    void get_Higgs_Nucleon_coupling_fN (Higgs_Nucleon_coupling_fN &result)
    {
      using namespace Pipes::get_Higgs_Nucleon_coupling_fN;

      const double sigmas = *Param["sigmas"]*1e-3, sigmal = *Param["sigmal"]*1e-3; // nuclear parameters in GeV (model input in MeV)
      const Spectrum SM = *Dep::SM_spectrum; // SM spectrum needed to get light quark masses

      const double z = 1.49; // isospin breaking ratio
      const double mu = SM.get(Par::mass1, "u_1"), md = SM.get(Par::mass1, "d_1"), ms = SM.get(Par::mass1, "d_2"); // light quark masses [GeV]
      const double mn = Gambit::m_neutron, mp = Gambit::m_proton; // nucleon masses [GeV]

      // intermediate quantities
      const double ml = 0.5*(mu+md);
      const double sigma0 = sigmal - sigmas*(2.*ml/ms);;
      const double y = 1 - sigma0/sigmal;

      std::vector<double> fu, fd, fs, mN = {mn, mp};

      for (size_t i(0); i<mN.size(); ++i)
      {
        fu.push_back(mu/(mu+md)*sigmal/mN[i]*(2*z+y*(1-z))/(1+z));
        fd.push_back(md/(mu+md)*sigmal/mN[i]*(2-y*(1-z))/(1+z));
        fs.push_back(ms/(mu+md)*sigmal/mN[i]*y);
      }

      result.neutron =  2./9. + 7./9.*(fu[0]+fd[0]+fs[0]);
      result.proton  =  2./9. + 7./9.*(fu[1]+fd[1]+fs[1]);
    }

    // Modified Inverse-Square Law (ISL) by adding a new Yukawa potential to the Newtonian gravitational potential: Vnew(r) = -(alpha*G*m1*m2)/r * exp(-r/lambda)
    // where alpha is the strenght of the new force and lambda its range

    // experimental parameters from Sushkov et al. 2011 arXiv:1108.2547
    const double rhoAu = 19, rhoTi = 4.5, rhog = 2.6, dAu = 700e-8, dTi = 100e-8, R = 15.6; // in cgs units

    // capability function returning the new force from the SuperRenormHP model for the experiment from Shuskov et al. 2011
    void New_Force_Sushkov2011_SuperRenormHP (daFunk::Funk &result)
    {
      using namespace Pipes::New_Force_Sushkov2011_SuperRenormHP;

      const double alpha = *Param["alpha"], lambda = *Param["lambda"];

      daFunk::Funk d = daFunk::var("d");

      daFunk::Funk force = 4*pow(pi, 2)*G_cgs*R*alpha*pow(lambda, 3)*exp(-d/lambda)*pow(rhoAu + (rhoTi-rhoAu)*exp(-dAu/lambda) + (rhog-rhoTi)*exp(-(dAu+dTi)/lambda), 2)*1e-5; // *1e-5 conversion from dyn(cgs) to N (SI)

      result = force;
    }

    // capability function to compute the likelihood from Sushkov et al. 2011
    void calc_lnL_ShortRangeForces_Sushkov2011 (double &result)
    {
      using namespace Pipes::calc_lnL_ShortRangeForces_Sushkov2011;

      daFunk::Funk ForceNew = *Dep::New_Force_Sushkov2011*1e12; // new force in pN

      static ASCIItableReader data = ASCIItableReader(GAMBIT_DIR "/DarkBit/data/ShortRangeForces/Sushkov2011.dat");
      data.setcolnames({"distance", "Fres", "sigma", "binWidth"});
      static std::vector<double> distance = data["distance"]; // [microns]
      static std::vector<double> Fres = data["Fres"]; // [pN]
      static std::vector<double> sigma = data["sigma"]; // [pN]
      static std::vector<double> width = data["binWidth"]; // [microns]

      std::vector<boost::shared_ptr<daFunk::FunkBase>> ForceNewBinned;
      std::vector<boost::shared_ptr<daFunk::FunkBound>> FnewBound;

      double d, delta;

      for (size_t i(0); i<distance.size(); ++i)
      {
        d = distance[i]*1e-4;
        delta = width[i]*1e-4;
        ForceNewBinned.push_back(ForceNew->gsl_integration("d", d-delta/2, d+delta/2)/delta);
        FnewBound.push_back(ForceNewBinned[i]->bind());
      }

      std::vector<double> likelihood;
      double norm, Fnew;

      for (size_t i(0); i<distance.size(); ++i)
      {
        norm = 1.; // we take the likelihood ratio to avoid having different normalizations accross the parameter space;
        Fnew = FnewBound[i]->eval();
        likelihood.push_back( (Fnew<Fres[i]) ? norm : norm*exp(-pow(Fres[i]-Fnew, 2)/pow(sigma[i], 2)) );
      }

      result = log(*std::min_element(likelihood.begin(), likelihood.end())); // we take the minimum likelihood, since we don't have the correlations between data bins
    }

  }
}
