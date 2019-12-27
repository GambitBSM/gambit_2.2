///  GAMBIT: Global and Modular BSM Inference Tool
///  *********************************************
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
#include <vector>
#include <string>
#include <iostream>

#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_monte_plain.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_min.h>

#include "gambit/Elements/gambit_module_headers.hpp"
#include "gambit/Utils/util_functions.hpp"
#include "gambit/Utils/numerical_constants.hpp"
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

    //------------- Physical constants and other useful things -------------// 

    const double pi=Gambit::pi;
    const double alpha=Gambit::alpha_EM; // fine structure constant
    const double v(246e9); // eV electroweak vev
    const double hbar=Gambit::hbar*1e9;  // eV.s
    const double cs=Gambit::s2cm*1e-2.; // cm/s speed of light
    const double s0(2891); // 1/cm³ current entropy density
    const double rho0(0.3e9); // eV/cm³ local DM density (Milky Way)
    const double rhoC(4.84e3); // ev/cm³ current critical density
    const double r0(26.2225e21); // cm Sun's distance from galactic center
    const double C(50./27.);
    const double sigma_gaussian(1e3);


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
      double m_H0;
      double m_t0;
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
    //         Support class to handle x-ray experiments              //
    ////////////////////////////////////////////////////////////////////

    //------------- Class declaration -------------//  

    class GammaRay
    {
      public:

      GammaRay(std::string experiment);
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
      ~GammaRay();

      protected:

      double m_J;
      double m_Emin;
      double m_Emax;
      double m_deltaOmega;
      double m_deltaE; // energy resolution
      std::vector<std::vector<double> > m_lRange;
      std::vector<std::vector<double> > m_bRange;
      std::string m_experiment;
      std::map<std::string, int> m_experimentMap;
    };

    // constructor
    GammaRay::GammaRay(std::string experiment) : m_experiment(experiment), m_experimentMap({{"INTEGRAL", 1}, {"HEAO", 2}})
    {
      switch(m_experimentMap[m_experiment])
      { 
        case 1 : 
          m_Emin = 20e3;
          m_Emax = 2e6;
          m_lRange = { {0., 60.} };
          m_bRange = { {75., 105.} };
          m_J = 3.65;
          m_deltaE = 8e3;
          break;

        case 2 :
          m_Emin = 4e3;
          m_Emax = 30e3;
          m_lRange = { {58., 109.}, {238., 289.} };
          m_bRange = { {0., 70.}, {110., 180.} }; 
          m_J = 3.88;
          m_deltaE = 4e3;
          break;

        default : 
          throw std::runtime_error("Wrong experiment name in GammaRay object");
          break;
      }
      set_deltaOmega();
    }

    // ----------- Functions to compute the solid angle of observation -----------
    
    // auxiliary function for gsl integration
    double deltaOmega (double x[], size_t dim, void *p)
    {
      (void)(p);
      (void)(dim);
      return sin(x[1]);
    }

    // computes the solide angle for a given galactic coordinates range (in degrees)
    double GammaRay::solidAngle(std::vector<double> lRange, std::vector<double> bRange)
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
    void GammaRay::set_deltaOmega()
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
    double GammaRay::flux(double const& E)
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
      GammaRay *experiment = static_cast<GammaRay*>(p);
      return experiment->flux(x);
    }

    // photon flux integrated over an interval deltaE, centered around E [photons/cm²/s]
    double GammaRay::fluxIntegrated(double const& E)
    {
      size_t neval; 
      double epsabs = 0.;
      double epsrel = 1e-2;
      double result, abserr;

      gsl_function F;
      F.function = &flux_gsl;
      F.params = this;

      gsl_integration_qng(&F, E-m_deltaE/2., E+m_deltaE/2., epsabs, epsrel, &result, &abserr, &neval);

      return result;
    }

    // standard deviation of the differential photon flux [photons/keV/cm²/s]
    double GammaRay::sigma(double const& E)
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
      GammaRay *experiment = static_cast<GammaRay*>(p);
      return pow(experiment->sigma(x), 2.);
    }

    // standard deviation of the photon flux integrated over an interval deltaE, centered around E [photons/cm²/s]
    double GammaRay::sigmaIntegrated(double const& E)
    {
      size_t neval; 
      double epsabs = 0.;
      double epsrel = 1e-2;
      double result, abserr;

      gsl_function F;
      F.function = &sigma_gsl;
      F.params = this;

      gsl_integration_qng(&F, E-m_deltaE/2., E+m_deltaE/2., epsabs, epsrel, &result, &abserr, &neval);

      return sqrt(result);
    }

    //------------- Elevator functions -------------// 

    double GammaRay::getDeltaOmega() const { return m_deltaOmega; }

    double GammaRay::getJ() const { return m_J; }

    double GammaRay::getEmin() const { return m_Emin; }

    double GammaRay::getEmax() const { return m_Emax; }

    double GammaRay::getDeltaE() const { return m_deltaE; }

    // destructor
    GammaRay::~GammaRay() { }

    ////////////////////////////////////////////////////////////////////
    //        Support class to handle the Higgs Portal model          //
    ////////////////////////////////////////////////////////////////////
    
    
    //------------- Class declaration -------------// 

    class HiggsPortal 
    {
      public:
      HiggsPortal(double const& theta, double const& ms, Cosmology const& cosmology);
      void setGamma();
      void setYs();
      void setParameters(double const& theta, double const& ms);
      double d2PhiEg(double const& E);
      double dPhiG(double const& E, double const& J);
      double gammaRayPrediction(double const& E, GammaRay const& experiment);
      double gammaRayPredictionIntegrated(double const& E, GammaRay const& experiment);
      double minimizeLikelihood(GammaRay& experiment);
      ~HiggsPortal();

      private:
      double m_theta;
      double m_ms;
      double m_gamma;
      double m_Ys;
      Cosmology m_cosmology;
    };

    // constructor
    HiggsPortal::HiggsPortal(double const& theta, double const& ms, Cosmology const& cosmology) : m_theta(theta), m_ms(ms), m_gamma(0.), m_Ys(0.), m_cosmology(cosmology) { setGamma(); setYs(); } 

    //------------- Functions to set the parameters of the model and related quantities -------------// 

    // decay rate to two photons [1/s]
    void HiggsPortal::setGamma() 
    {
      m_gamma = (m_theta*m_theta*alpha*alpha*m_ms*m_ms*m_ms*C*C)/(256.*pi*pi*pi*v*v)/hbar;  
    }

    // initial abundance [K/eV]
    void HiggsPortal::setYs()
    {
      m_Ys = 4e11*m_theta*m_theta;
    }

    // sets new values for the parameters of the model and updates related quantities
    void HiggsPortal::setParameters(double const& theta, double const& ms)
    {
      m_theta = theta;
      m_ms = ms;
      setGamma();
      setYs();
    }

    //------------- Functions to compute the predicted photon flux -------------// 

    // useful structure
    struct gammaRayLikelihood_params {HiggsPortal *model; GammaRay experiment;};

    // cosmological contribution to the differential photon flux [photons/eV/cm²/s/sr]
    double HiggsPortal::d2PhiEg(double const& E) 
    {
      double x = m_ms/2./E;
      double z = x - 1.;

      return 2.*1./(4*pi)*(m_gamma*m_Ys*s0*m_ms*cs*exp(-m_gamma*m_cosmology.age(z)[0]))/(m_ms*m_cosmology.getH0()*E)/sqrt( m_cosmology.getOmegaM()*pow(x, 3.) + m_cosmology.getOmegaLambda() + m_cosmology.getOmegaK()*pow(x, 2.) + m_cosmology.getOmegaR()*pow(x, 4.) );
    }

    // galactic (Milky Way) contribution to the differential photon flux [photons/eV/cm²/s]
    double HiggsPortal::dPhiG(double const& E, double const& J)
    {
      return 2.*(r0*rho0*m_gamma*J*m_Ys*s0*m_ms*exp(-m_cosmology.get_t0()*m_gamma))/(4.*pi*m_ms*m_cosmology.getOmegaDM()*rhoC)/sqrt(2*pi*sigma_gaussian*sigma_gaussian)*exp(-pow(E-m_ms/2.,2)/(2*sigma_gaussian*sigma_gaussian));
    }

    // total predicted differential photon flux for a given x-ray experiment [photons/eV/cm²/s]
    double HiggsPortal::gammaRayPrediction(double const& E, GammaRay const& experiment)
    {
      double J, deltaOmega;

      J = experiment.getJ();
      deltaOmega = experiment.getDeltaOmega();
      return dPhiG(E, J) + (d2PhiEg(E) * deltaOmega);
    }

    // auxiliary function for gsl integration
    double gammaRayPredictionIntegrated_gsl(double x, void *p)
    {
      gammaRayLikelihood_params *params = static_cast<gammaRayLikelihood_params*>(p);
      HiggsPortal *model = params->model;
      GammaRay experiment = params->experiment;
  
      return model->gammaRayPrediction(x, experiment);
    }

    // total predicted photon flux integrated over an interval deltaE, centered around E [photons/cm²/s]
    double HiggsPortal::gammaRayPredictionIntegrated(double const& E, GammaRay const& experiment)
    {
      size_t n = 1e4; 

      gsl_integration_workspace *w =  gsl_integration_workspace_alloc(n);

      double epsabs = 0.;
      double epsrel = 1e-2;
      size_t limit = 1e3;
      double result, abserr;
      int key = 6;
      double deltaE = experiment.getDeltaE();

      gammaRayLikelihood_params params = {this, experiment};

      gsl_function F;
      F.function = &gammaRayPredictionIntegrated_gsl;
      F.params = &params;

      gsl_integration_qag(&F, E-deltaE/2., E+deltaE/2., epsabs, epsrel, limit, key, w, &result, &abserr);

      gsl_integration_workspace_free(w);

      return result;
    }

    //------------- Functions to compute x-ray likelihoods -------------// 

    // auxiliary function for gsl minimization
    double gammaRayLikelihood(double E, void * p) 
    {
      gammaRayLikelihood_params *params = static_cast<gammaRayLikelihood_params*>(p);
      HiggsPortal *model = params->model;
      GammaRay experiment = params->experiment;

      double data = experiment.fluxIntegrated(E);
      double sigma = experiment.sigmaIntegrated(E);
      double prediction = model->gammaRayPredictionIntegrated(E, experiment);

      if (prediction>=data) 
      {  
        return exp(-pow(data-prediction,2.)/(2.*sigma*sigma));
      }

      else 
      {
        return 1.;
      }
    }

    // computes the energy E which minimizes the likelihood
    double HiggsPortal::minimizeLikelihood(GammaRay& experiment) {

      int status;
      int iter = 0, max_iter = 100;
      const gsl_min_fminimizer_type *T;
      gsl_min_fminimizer *s;
      double a = experiment.getEmin(), b = m_ms/2.;
      double m = (a+b)/2.;
      gsl_function F;

      gammaRayLikelihood_params params = {this, experiment};

      F.function = &gammaRayLikelihood;
      F.params = &params;
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

    ////////////////////////////////////////////////////////////////////
    //                   Capability functions                         //
    ////////////////////////////////////////////////////////////////////

    // instance of the cosmology class
    auto cosmology = Cosmology();
    const double t0 = cosmology.get_t0();

    // capability function to provide the decay rate to two photons    
    void calc_SuperRenormalizableHiggsPortalDM_decay_rate(double &result)
    {
      using namespace Pipes::calc_SuperRenormalizableHiggsPortalDM_decay_rate;
      double mS = *Param["mS"], theta = *Param["theta"];
      
      result = (theta*theta*alpha*alpha*mS*mS*mS*C*C)/(256.*pi*pi*pi*v*v)/hbar;
    }
    
    void calc_SuperRenormalizableHiggsPortal

  }
}
    
      
    

