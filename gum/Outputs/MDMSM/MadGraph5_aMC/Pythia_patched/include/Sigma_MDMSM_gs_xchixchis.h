//==========================================================================
// This file has been automatically generated for Pythia 8
// MadGraph5_aMC@NLO v. 2.8.3.2, 2021-02-02
// By the MadGraph5_aMC@NLO Development Team
// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
//==========================================================================

#ifndef Pythia8_Sigma_MDMSM_gs_xchixchis_H
#define Pythia8_Sigma_MDMSM_gs_xchixchis_H

#include <complex> 

#include "Pythia8/SigmaProcess.h"
#include "Parameters_MDMSM.h"

using namespace std; 

namespace Pythia8 
{
//==========================================================================
// A class for calculating the matrix elements for
// Process: g s > ~chi ~chi s NP<=2 WEIGHTED<=3 @2
// Process: g s~ > ~chi ~chi s~ NP<=2 WEIGHTED<=3 @2
//--------------------------------------------------------------------------

class Sigma_MDMSM_gs_xchixchis : public Sigma3Process 
{
  public:

    // Constructor.
    Sigma_MDMSM_gs_xchixchis() {}

    // Initialize process.
    virtual void initProc(); 

    // Calculate flavour-independent parts of cross section.
    virtual void sigmaKin(); 

    // Evaluate sigmaHat(sHat).
    virtual double sigmaHat(); 

    // Select flavour, colour and anticolour.
    virtual void setIdColAcol(); 

    // Evaluate weight for decay angles.
    virtual double weightDecay(Event& process, int iResBeg, int iResEnd); 

    // Info on the subprocess.
    virtual string name() const {return "gs_xchixchis (MDMSM)";}

    virtual int code() const {return 10206;}

    virtual string inFlux() const {return "qg";}
    int id3Mass() const {return 52;}
    int id4Mass() const {return 52;}
    int id5Mass() const {return 3;}
    virtual int resonanceA() const {return 99902;}
    // Tell Pythia that sigmaHat returns the ME^2
    virtual bool convertM2() const {return true;}

  private:

    // Private functions to calculate the matrix element for all subprocesses
    // Calculate wavefunctions
    void calculate_wavefunctions(const int perm[], const int hel[]); 
    static const int nwavefuncs = 10; 
    std::complex<double> w[nwavefuncs][18]; 
    static const int namplitudes = 2; 
    std::complex<double> amp[namplitudes]; 
    double matrix_2_gs_xchixchis(); 
    double matrix_2_gsx_xchixchisx(); 

    // Constants for array limits
    static const int nexternal = 5; 
    static const int nprocesses = 4; 

    // Store the matrix element value from sigmaKin
    double matrix_element[nprocesses]; 

    // Color flows, used when selecting color
    double * jamp2[nprocesses]; 

    // Pointer to the model parameters
    Parameters_MDMSM * pars; 

}; 

}  // end namespace Pythia8

#endif  // Pythia8_Sigma_MDMSM_gs_xchixchis_H

