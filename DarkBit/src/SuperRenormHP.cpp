///  \file
///
///  Super Renormalizable Higgs Portal DM specific module functions for DarkBit
///
///  *********************************************
///
///  Authors (add name and date if you modify):
///
///  \author Iñigo Saez Casares
///  \date 2019 December
///
///  *********************************************

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

#include "fjcore.hh"
#include "gambit/DecayBit/DecayBit_rollcall.hpp"
#include "gambit/Elements/gambit_module_headers.hpp"
#include "gambit/Utils/util_functions.hpp"
#include "gambit/Utils/ascii_table_reader.hpp"
#include "gambit/Utils/numerical_constants.hpp"
#include "gambit/Utils/statistics.hpp"
#include "gambit/DarkBit/DarkBit_rollcall.hpp"

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
    const double alphaS = pi; // strong coupling constant
    const double v = 246e9; // electroweak vev [eV]
    const double hbar = Gambit::hbar*1e9;  // reduced Planck constant [eV.s]
    const double cs = Gambit::s2cm; // speed of light [cm/s]
    const double kb = Gambit::K2eV; // Boltzmann constant [eV/K]
    const double G = 6.674e-8; // Gravitational constant [cm³/g/s²]
    const double C(50./27.); // loop function from the decay of the Higgs boson into two photons

    // cosmological constants
    const double s0(2891); // current entropy density [1/cm³]
    const double rhoC(4.84e3); // current critical density [eV/cm³]

    // astrophysical constants
    const double r0(26.2225e21); // Sun's distance from the galactic center [cm]
    const double Rsun(5.9598e10); // Solar radius [cm]
    const double L0 = 2.388672e45; // solar photon luminosity [eV/s]
    const double rho0(0.3e9); // local DM density (Milky Way) [eV/cm³]

    // Minimum finite result returnable from log(double x);
    const double logmin = log(std::numeric_limits<double>::min());

    ////////////////////////////////////////////////////////////////////
    //      Support class to compute cosmological observables         //
    ////////////////////////////////////////////////////////////////////


    //------------- Class declaration -------------// 

    class Cosmology 
    {
      public :

        Cosmology();

        void set_t0();

        double getOmegaK() const;
        double getOmegaDM() const;
        double getOmegaM() const;
        double getOmegaLambda() const;
        double getOmegaR() const;
        double getOmegaB() const;
        double getH0() const;
        double get_t0() const;

        std::vector<double> age (double redshift);

        ~Cosmology();

      protected :

        double m_OmegaK;
        double m_OmegaDM;
        double m_OmegaM;
        double m_OmegaLambda;
        double m_OmegaR;
        double m_OmegaB;
        double m_H0; // Hubble constant
        double m_t0; // age of the universe at redshift zero
    };

    // constructor
    Cosmology::Cosmology() : m_OmegaK(-0.005), m_OmegaDM(0.285), m_OmegaM(0.308), m_OmegaLambda(0.692), m_OmegaR(0), m_OmegaB(0.308-0.285), m_H0(67.81e-19/3.085), m_t0(0) { set_t0(); }

    //------------- Functions to compute the age of the Universe at a given redshift -------------// 

    // auxiliary function for gsl integration
    double age_f (double x, void *p)
    {
      Cosmology *cosmology = static_cast<Cosmology*>(p);
      double OmegaK = cosmology->getOmegaK();
      double OmegaM = cosmology->getOmegaM();
      double OmegaLambda = cosmology->getOmegaLambda();
      double OmegaR = cosmology->getOmegaR();

      return 1./sqrt( OmegaM*pow(1+x, 5.) + OmegaLambda*pow(1+x, 2.) + OmegaK*pow(1+x, 4.) + OmegaR*pow(1+x, 6.) );
    }

    // computes the age of the Universe at a given redshift ([0] age, [1] abserr)
    std::vector<double> Cosmology::age (double redshift)
    {
      size_t n = 1e4; 

      gsl_integration_workspace *w =  gsl_integration_workspace_alloc(n);

      double epsabs = 0.;
      double epsrel = 1e-3;
      size_t limit = 1e3;
      double result, abserr;

      gsl_function F;
      F.function = &age_f;
      F.params = this;

      gsl_integration_qagiu(&F, redshift, epsabs, epsrel, limit, w, &result, &abserr);

      gsl_integration_workspace_free(w);

      return {result/m_H0, abserr/m_H0};
    }

    //------------- Elevator functions -------------// 

    void Cosmology::set_t0 () { m_t0 = age(0)[0]; }
    double Cosmology::getOmegaK () const { return m_OmegaK; }
    double Cosmology::getOmegaDM () const { return m_OmegaDM; }
    double Cosmology::getOmegaM () const { return m_OmegaM; }
    double Cosmology::getOmegaLambda () const { return m_OmegaLambda; }
    double Cosmology::getOmegaR () const { return m_OmegaR; }
    double Cosmology::getOmegaB () const { return m_OmegaB; }
    double Cosmology::getH0 () const { return m_H0; }
    double Cosmology::get_t0 () const { return m_t0; }

    // destructor
    Cosmology::~Cosmology() {}

    ////////////////////////////////////////////////////////////////////
    //         Support class to handle X-ray experiments              //
    ////////////////////////////////////////////////////////////////////

    //------------- Class declaration -------------//  

    class Xray
    {
      public:

        Xray(std::string experiment);
        double solidAngle(std::vector<double> lRange, std::vector<double> bRange);
        void set_deltaOmega();
        double getDeltaOmega() const;
        double getJ() const;
        double getEmin() const;
        double getEmax() const;
        double getDeltaE() const;
        double flux(double const& E);
        double sigma(double const& E);
        double fluxIntegrated(double const& E);
        double sigmaIntegrated(double const& E);
        double deltaE(double const& E) const;
        ~Xray();

      protected:

        double m_J;
        double m_Emin;
        double m_Emax;
        double m_deltaOmega;
        double m_deltaE; // energy resolution in percentage of the energy 
        std::vector<std::vector<double> > m_lRange;
        std::vector<std::vector<double> > m_bRange;
        std::string m_experiment;
        std::map<std::string, int> m_experimentMap;
    };

    // constructor
    Xray::Xray(std::string experiment) : m_experiment(experiment), m_experimentMap({{"INTEGRAL", 1}, {"HEAO", 2}})
    {
      switch(m_experimentMap[m_experiment])
      { 
        case 1 : 
          m_Emin = 20e3;
          m_Emax = 2e6;
          m_lRange = { {0., 60.} };
          m_bRange = { {75., 105.} };
          m_J = 3.65;
          m_deltaE = 0.1;
          m_deltaOmega = 0.542068;
          break;

        case 2 :
          m_Emin = 4e3;
          m_Emax = 30e3;
          m_lRange = { {58., 109.}, {238., 289.} };
          m_bRange = { {0., 70.}, {110., 180.} }; 
          m_J = 3.88;
          m_deltaE = 0.1;
          m_deltaOmega = 1.17135;
          break;

        default : 
          throw std::runtime_error("Wrong experiment name in Xray object");
          break;
      }
      //set_deltaOmega();
    }

    //------------- Function returning the energy dispersion of the instrument -------------// 

    double Xray::deltaE (double const& E) const
    {
      return m_deltaE*E;
    }

    // ----------- Functions to compute the solid angle of observation -------------// 

    // auxiliary function for gsl integration
    double deltaOmega (double x[], size_t dim, void *p)
    {
      (void)(p);
      (void)(dim);
      return sin(x[1]);
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
          return 4.8e-8*pow(E/100e3,-1.55) + 6.6e-8*exp(-(E-50e3)/7.5e3);
          break;

        case 2 :
          return 7.877*pow(10.,0.87)*pow(E,-0.29)*exp(-(E/41.13e3))/E*m_deltaOmega;
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

    const double int_factor(1.1); // integration width = int_factor*energy dispersion instrument

    // photon flux integrated over an interval deltaE, centered around E [photons/cm²/s]
    double Xray::fluxIntegrated(double const& E)
    {
      size_t neval; 
      double epsabs = 0.;
      double epsrel = 1e-2;
      double result, abserr;
      double delta = int_factor*deltaE(E);

      gsl_function F;
      F.function = &flux_gsl;
      F.params = this;

      gsl_integration_qng(&F, E-delta/2., E+delta/2., epsabs, epsrel, &result, &abserr, &neval);

      return result;
    }

    // standard deviation of the differential photon flux [photons/keV/cm²/s]
    double Xray::sigma(double const& E)
    {

      switch(m_experimentMap[m_experiment]) 
      {
        case 1 : 
          return sqrt(pow(pow(E/100e3,-1.55),2)*pow(0.6e-5,2) + pow(4.8e-8*1.55*pow(E/100e3,-2.55),2)*pow(0.25,2) + exp(-2*(E-50e3)/7.5e3)*pow(0.5e-8, 2.) + pow(6.6e-8, 2.)*pow((E-50e3)/pow(7.5e3, 2.), 2.)*exp(-2*(E-50e3)/7.5e3)*pow(1e3, 2.)); 
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

    // standard deviation of the photon flux integrated over an interval deltaE, centered around E [photons/cm²/s]
    double Xray::sigmaIntegrated(double const& E)
    {
      size_t neval; 
      double epsabs = 0.;
      double epsrel = 1e-2;
      double result, abserr;
      double delta = int_factor*deltaE(E);

      gsl_function F;
      F.function = &sigma_gsl;
      F.params = this;

      gsl_integration_qng(&F, E-delta/2., E+delta/2., epsabs, epsrel, &result, &abserr, &neval);

      return sqrt(result);
    }

    //------------- Elevator functions -------------// 

    double Xray::getDeltaOmega() const { return m_deltaOmega; }

    double Xray::getJ() const { return m_J; }

    double Xray::getEmin() const { return m_Emin; }

    double Xray::getEmax() const { return m_Emax; }

    double Xray::getDeltaE() const { return m_deltaE; }

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
        wp.push_back(sqrt(4*pi*alphaEM*sumz*cs*cs*cs*hbar*hbar*hbar/me));
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
        result += D/N*F*pow(cs, 6)*pow(hbar, 6);
      }

      // uses a first order development at high x instead of the full functions
      catch (int gsl_errno)
      {
        double F = sqrt(pi/2./x)*(1-(1/8./x));
        result += D/N*F*pow(cs, 6)*pow(hbar, 6);
      }

      if (w > wp) { result += 8*pi*pow(alphaEM, 2)*ne*sqrt(1-pow(wp/w, 2))/pow(me, 3)/3.*pow(cs, 3)*pow(hbar, 3); }

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

      return 4*pi*pow(Rsun, 3)*(result1+result2)/pow(cs, 3)/pow(hbar, 4);
    }

    void StellarModel::Ls_interpolate ()
    {
      const int nPoints = 150;
      const double mMin = 5e-1, mMax = mSmax;
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


    //------------- Functions to compute X-ray likelihoods -------------// 

    // useful structure
    struct XrayLikelihood_params {double mass; double gamma; double density; Xray experiment; Cosmology cosmology;};

    // cosmological contribution to the differential photon flux [photons/eV/cm²/s/sr]
    double dPhiEg(double const& E, XrayLikelihood_params *params) 
    {
      double mass = params->mass, gamma = params->gamma, density = params->density;
      Xray experiment = params->experiment;
      Cosmology cosmology = params->cosmology;

      double x = mass/2./E;
      double z = x - 1.;

      return experiment.getDeltaOmega()*2.*1./(4*pi)*(gamma*density*cs*exp(-gamma*cosmology.age(z)[0]))/(mass*cosmology.getH0()*E)/sqrt( cosmology.getOmegaM()*pow(x, 3.) + cosmology.getOmegaLambda() + cosmology.getOmegaK()*pow(x, 2.) + cosmology.getOmegaR()*pow(x, 4.) );
    }


    // galactic (Milky Way) contribution to the differential photon flux [photons/eV/cm²/s]
    const double s(0.5); // sigma of the gaussian for the galactic emission line = s*energy dispersion instrument

    double dPhiG(double const& E, XrayLikelihood_params *params)
    {
      double mass = params->mass, gamma = params->gamma, density = params->density;
      Xray experiment = params->experiment;
      Cosmology cosmology = params->cosmology;

      double sigma = s*experiment.deltaE(E); // standard deviation of the gaussian modelling the enery dispersion of the instrument
      return 2.*(r0*rho0*gamma*experiment.getJ()*density*exp(-cosmology.get_t0()*gamma))/(4.*pi*mass*cosmology.getOmegaDM()*rhoC)/sqrt(2*pi*sigma*sigma)*exp(-pow(E-mass/2.,2)/(2*sigma*sigma));
    }

    // total predicted differential photon flux for a given X-ray experiment [photons/eV/cm²/s]
    double XrayPrediction(double const& E, XrayLikelihood_params *params)
    {
      return dPhiG(E, params) + (dPhiEg(E, params));
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

    // auxiliary function for gsl minimization
    double XrayLikelihood(double E, void * p)
    {
      XrayLikelihood_params *params = static_cast<XrayLikelihood_params*>(p);
      Xray experiment = params->experiment;

      double data = experiment.fluxIntegrated(E);
      double sigma = experiment.sigmaIntegrated(E);
      double prediction = XrayPredictionIntegrated(E, params);

      return (prediction>=data) ? exp(-pow(data-prediction,2.)/(2.*sigma*sigma)) : 1.;
    }

    // computes the energy E which minimizes the likelihood
    double minimizeLikelihood(XrayLikelihood_params *params)
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

      F.function = &XrayLikelihood;
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

    // capability function to provide the initial energy density as produced by freeze-in
    void SuperRenormHP_initial_density (double &result)
    {
      using namespace Pipes::SuperRenormHP_initial_density;
      double mS = *Param["mS"], theta = *Param["theta"];
      result = 4e11*theta*theta*s0*mS;
    }

    // capability function to provide the decay rate to two photons
    void SuperRenormHP_decay_rate_2photons (double &result)
    {
      using namespace Pipes::SuperRenormHP_decay_rate_2photons;
      double mS = *Param["mS"], theta = *Param["theta"];

      result = (theta*theta*alphaEM*alphaEM*mS*mS*mS*C*C)/(256.*pi*pi*pi*v*v)/hbar;
    }

    // capability function to provide the mass of the DM candidate
    void SuperRenormHP_mass (double &result)
    {
      using namespace Pipes::SuperRenormHP_mass;

      result = *Param["mS"];
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

    // capability function to compute the X-ray Likelihood from the INTEGRAL experiment
    void calc_lnL_INTEGRAL(double &result)
    {
      using namespace Pipes::calc_lnL_INTEGRAL;

      static Cosmology cosmology = Cosmology();

      static Xray experiment = Xray("INTEGRAL");

      double mass = *Dep::DM_mass, gamma = *Dep::decay_rate_2photons, density = *Dep::initial_density;

      XrayLikelihood_params params = {mass, gamma, density, experiment, cosmology};

      double Emin = experiment.getEmin(), Emax = experiment.getEmax(), E, lik1, lik2;

      if (mass >= 1e6) { result = 0; }

      else if (mass >= 2.*Emin)
      {
        // modifies the gsl error handler and stores the default one
        gsl_error_handler_t *old_handler = gsl_set_error_handler (&handler);
        try
        {
          E = minimizeLikelihood(&params);
          result = log(XrayLikelihood(E, &params));
        }

        catch (int gsl_errno)
        {
          lik1 = XrayLikelihood(Emin+experiment.deltaE(Emin), &params);
          lik2 = XrayLikelihood(fmin(mass/2., Emax-experiment.deltaE(Emax)), &params);
          result = log(fmin(lik1, lik2));
        }
        // restores the default gsl error handler
        gsl_set_error_handler (old_handler);
      }

      else { result = 0; }
    }

    // capability function to compute the X-ray Likelihood from the HEAO-1 A2 experiment
    void calc_lnL_HEAO(double &result)
    {
      using namespace Pipes::calc_lnL_HEAO;

      static Cosmology cosmology = Cosmology();

      static Xray experiment = Xray("HEAO");

      const double mass = *Dep::DM_mass, gamma = *Dep::decay_rate_2photons, density = *Dep::initial_density;

      XrayLikelihood_params params = {mass, gamma, density, experiment, cosmology};

      const double Emin = experiment.getEmin(), Emax = experiment.getEmax();
      double E, lik1, lik2;

      if (mass >= 1e6) { result = 0; }

      else if (mass >= 2.*Emin)
      {
        // modifies the gsl error handler and stores the default one
        gsl_error_handler_t *old_handler = gsl_set_error_handler (&handler);
        try
        {
          E = minimizeLikelihood(&params);
          result = log(XrayLikelihood(E, &params));
        }

        catch (int gsl_errno)
        {
          lik1 = XrayLikelihood(Emin+experiment.deltaE(Emin), &params);
          lik2 = XrayLikelihood(fmin(mass/2., Emax-experiment.deltaE(Emax)), &params);
          result = log(fmin(lik1, lik2));
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
      const double mS = *Param["mS"], theta = *Param["theta"];

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
    //
    // capability function to compute the predicted solar neutrino flux (Be7)
    void SuperRenormHP_solar_neutrino_flux_Be7 (double &result)
    {
      using namespace Pipes::SuperRenormHP_solar_neutrino_flux_Be7;

      const double Ls = *Dep::solar_DM_luminosity;
      const double alpha = runOptions->getValueOrDef<double>(4., "alpha");

      const double Phi0 = 4.71e9;

      result = Phi0*pow(1+Ls/L0, alpha);
    }

    // capability function to compute the likelihood from solar B8 neutrino flux
    void calc_lnL_solar_neutrino_B8 (double &result)
    { 
      using namespace Pipes::calc_lnL_solar_neutrino_B8;

      const bool profile = runOptions->getValueOrDef<bool>(false, "profile_systematics");
      const double Phi_predicted = *Dep::solar_neutrino_flux_B8;

      const double Phi_obs = 5e6;
      const double sigma_obs = 0.03*Phi_obs, sigma_theo = 0.14*Phi_predicted;

      result = Stats::gaussian_upper_limit(Phi_predicted, Phi_obs, sigma_theo, sigma_obs, profile);
    }

    // capability function to compute the likelihood from solar Be7 neutrino flux
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

    // Modified Inverse-Square Law (ISL) by adding a new Yukawa potential to the Newtonian gravitational potential: Vnew(r) = -(alpha*G*m1*m2)/r * exp(-r/lambda)
    // where alpha is the strenght of the new force and lambda its range

    // capability function returning the values of alpha and lambda for the SuperRenormHP model
    void SuperRenormHP_ISL_Yukawa (std::vector<double> &result) 
    {
      using namespace Pipes::SuperRenormHP_ISL_Yukawa;

      const double mS = *Param["mS"], theta = *Param["theta"];
      const double tau = 4*pow(mT/mH, 2.);
      const double g = alphaS/8./pi*tau*(1+(1-tau)*pow(asin(sqrt(1./tau)), 2.)); // effective Higgs-gluon coupling constant (generated by top quark loop)

      const double f = v/theta/g;
      result = {4*pi*pow(Mp/f, 2), hbar*cs/mS}; // lambda is given in cm
    }

    // experimental parameters from Sushkov et al. 2011 arXiv:1108.2547
    const double rhoAu = 19, rhoTi = 4.5, rhog = 2.6, dAu = 700e-8, dTi = 100e-8, R = 15.6; // in cgs units

    // capability function to compute the likelihood from Sushkov et al. 2011
    void calc_lnL_ShortRangeForces_Sushkov2011 (double &result)
    {
      using namespace Pipes::calc_lnL_ShortRangeForces_Sushkov2011;

      const std::vector<double> parameters = *Dep::ISL_Yukawa;
      const double alpha = parameters[0], lambda = parameters[1];

      ASCIItableReader data = ASCIItableReader(GAMBIT_DIR "/DarkBit/data/ShortRangeForces/Sushkov2011.dat");
      data.setcolnames({"distance", "Fres", "sigma", "binWidth"});
      std::vector<double> distance = data["distance"]; // [microns]
      std::vector<double> Fres = data["Fres"]; // [pN]
      std::vector<double> sigma = data["sigma"]; // [pN]
      std::vector<double> width = data["binWidth"]; // [microns]

      std::vector<double> Fnew;
      const double factor = 4*pi*G*R*alpha*pow(lambda, 3)*lambda*pow(rhoAu + (rhoTi-rhoAu)*exp(-dAu/lambda) + (rhog-rhoTi)*exp(-(dAu+dTi)/lambda), 2)*1e-5*1e12;

      double d, delta;

      for (size_t i(0); i<distance.size(); ++i)
      {
        d = distance[i]*1e-4;
        delta = width[i]*1e-4;
        Fnew.push_back(factor/delta*(exp(-(d-delta/2.)/lambda)-exp(-(d+delta/2.)/lambda))); // new force in pN
      }

      std::vector<double> likelihood;
      double norm;

      for (size_t i(0); i<distance.size(); ++i)
      {
        norm = 1./sqrt(2*pi)/sigma[i];
        likelihood.push_back( (Fnew[i]<Fres[i]) ? norm : norm*exp(-pow(Fres[i]-Fnew[i], 2)/pow(sigma[i], 2)) );
      }

      result = log(*std::min_element(likelihood.begin(), likelihood.end()));
    }

  }
}
