//==========================================================================
// This file has been automatically generated for Pythia 8
// MadGraph5_aMC@NLO v. 2.8.3.2, 2021-02-02
// By the MadGraph5_aMC@NLO Development Team
// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
//==========================================================================

#ifndef Pythia8_parameters_MDMSM_H
#define Pythia8_parameters_MDMSM_H

#include <complex> 

#include "Pythia8/ParticleData.h"
#include "Pythia8/StandardModel.h"
#include "Pythia8/SusyLesHouches.h"

using namespace std; 

using namespace Pythia8; 

class Parameters_MDMSM
{
  public:

    static Parameters_MDMSM * getInstance(); 

    // Model parameters independent of aS
    std::complex<double> mdl_complexi, mdl_I1b11, mdl_I1b22, mdl_I1b33,
        mdl_I2b11, mdl_I2b22, mdl_I2b33, mdl_I3b11, mdl_I3b22, mdl_I3b33,
        mdl_I4b11, mdl_I4b22, mdl_I4b33;
    double mdl_wY, mdl_Wchi, mdl_WH, mdl_WT, mdl_WW, mdl_WZ, mdl_mY, mdl_mchi,
        mdl_MH, mdl_MB, mdl_MS, mdl_MD, mdl_MT, mdl_MC, mdl_MU, mdl_MTA,
        mdl_MMU, mdl_Me, mdl_MZ, mdl_ymtau, mdl_ymm, mdl_yme, mdl_ymt, mdl_ymb,
        mdl_ymc, mdl_yms, mdl_ymup, mdl_ymdo, mdl_Gf, aEWM1, mdl_gggY2,
        mdl_gggY1, mdl_cY, mdl_gchi, ZERO, mdl_MZ__exp__2, mdl_MZ__exp__4,
        mdl_sqrt__2, mdl_MH__exp__2, mdl_aEW, mdl_MW, mdl_sqrt__aEW, mdl_ee,
        mdl_MW__exp__2, mdl_sw2, mdl_cw, mdl_sqrt__sw2, mdl_sw, mdl_g1, mdl_gw,
        mdl_vev, mdl_vev__exp__2, mdl_lam, mdl_yb, mdl_yc, mdl_ydo, mdl_ye,
        mdl_ym, mdl_ys, mdl_yt, mdl_ytau, mdl_yup, mdl_muH, mdl_ee__exp__2,
        mdl_sw__exp__2, mdl_cw__exp__2;
    // Model parameters dependent on aS
    double aS, mdl_sqrt__aS, G, mdl_G__exp__2; 
    // Model couplings independent of aS
    std::complex<double> GC_13, GC_14, GC_1, GC_2, GC_3, GC_4, GC_5, GC_6,
        GC_7, GC_8, GC_9, GC_17, GC_18, GC_19, GC_20, GC_21, GC_22, GC_23,
        GC_24, GC_25, GC_26, GC_27, GC_28, GC_29, GC_30, GC_31, GC_32, GC_33,
        GC_34, GC_35, GC_36, GC_37, GC_38, GC_39, GC_40, GC_41, GC_42, GC_43,
        GC_44, GC_45, GC_46, GC_47, GC_48, GC_49, GC_50, GC_51, GC_52, GC_53,
        GC_54, GC_55, GC_56, GC_57, GC_58, GC_59, GC_60, GC_61, GC_62, GC_63,
        GC_64, GC_65, GC_66, GC_67, GC_68, GC_69, GC_70, GC_71, GC_72, GC_73,
        GC_74, GC_75, GC_76, GC_77, GC_78, GC_79, GC_80, GC_81, GC_82, GC_83,
        GC_84, GC_85, GC_86, GC_87, GC_88, GC_89, GC_90, GC_91, GC_92, GC_93,
        GC_94, GC_95, GC_96, GC_97, GC_98, GC_99, GC_100, GC_101, GC_102;
    // Model couplings dependent on aS
    std::complex<double> GC_15, GC_11, GC_12, GC_16, GC_10; 

    // Set parameters that are unchanged during the run
    void setIndependentParameters(ParticleData * & pd, Couplings * & csm,
        SusyLesHouches * & slhaPtr);
    // Set couplings that are unchanged during the run
    void setIndependentCouplings(); 
    // Set parameters that are changed event by event
    void setDependentParameters(ParticleData * & pd, Couplings * & csm,
        SusyLesHouches * & slhaPtr, double alpS);
    // TMP: hardcoded bogus implementation with no arguments since this
    // is being called from within the matrix elements.
    void setDependentParameters() {}; 

    // Set couplings that are changed event by event
    void setDependentCouplings(); 

    // Print parameters that are unchanged during the run
    void printIndependentParameters(); 
    // Print couplings that are unchanged during the run
    void printIndependentCouplings(); 
    // Print parameters that are changed event by event
    void printDependentParameters(); 
    // Print couplings that are changed event by event
    void printDependentCouplings(); 

  private:
    static Parameters_MDMSM * instance; 
}; 

#endif  // Pythia8_parameters_MDMSM_H

