#ifndef __enum_decl_copies_FlexibleSUSY_CMSSM_2_0_1_h__
#define __enum_decl_copies_FlexibleSUSY_CMSSM_2_0_1_h__


#include "identification.hpp"

namespace CAT_3(BACKENDNAME,_,SAFE_VERSION)
{
    
    
    namespace flexiblesusy
    {
        typedef enum {kVerbose=0, kDebug=1, kInfo=2, kWarning=3, kError=4, kFatal=5} ELogLevel;
    }
    namespace softsusy
    {
        typedef enum {mUp=1, mCharm=2, mTop=3, mDown=4, mStrange=5, mBottom=6, mElectron=7, mMuon=8, mTau=9} mass;
    }
    namespace softsusy
    {
        typedef enum {ALPHA=1, ALPHAS=2} leGauge;
    }
    namespace softsusy
    {
        typedef enum {alpha_em_MSbar_at_MZ=0, alpha_s_MSbar_at_MZ=1, GFermi=2, MZ_pole=3, MW_pole=4, Mv1_pole=5, Mv2_pole=6, Mv3_pole=7, Me_pole=8, Mm_pole=9, Mtau_pole=10, mu_2GeV=11, ms_2GeV=12, Mt_pole=13, md_2GeV=14, mc_mc=15, mb_mb=16, CKM_theta_12=17, CKM_theta_13=18, CKM_theta_23=19, CKM_delta=20, PMNS_theta_12=21, PMNS_theta_13=22, PMNS_theta_23=23, PMNS_delta=24, PMNS_alpha_1=25, PMNS_alpha_2=26, NUMBER_OF_LOW_ENERGY_INPUT_PARAMETERS=27} QedQcd_input_parmeters;
    }
    
    
}


#include "gambit/Backends/backend_undefs.hpp"

#endif /* __enum_decl_copies_FlexibleSUSY_CMSSM_2_0_1_h__ */
