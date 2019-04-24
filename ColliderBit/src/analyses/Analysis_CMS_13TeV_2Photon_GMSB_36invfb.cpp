///
///  \author Anders Kvellestad
///  \date 2019 Apr
///
///  *********************************************

/* 
  Based on: 
    "Search for supersymmetry in final states with photons and missing transverse momentum in proton-proton collisions at 13 TeV" 

    http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-17-011/index.html
    arxiv:1903.07070

  Notes:


  Event selection:

    Apply photon detector efficiency?

    Trigger:
    - Two photons
    - Leading photon: pT > 30 
    - Subleading photon: pT > 18
    - Invariant mass m_gg > 95 
    - Isolation and cluster shape requirements (not implemented)

    Photon selection:
    - pT > 40
    - |eta| < 1.44
    - Isolated from other reconstructed particles by 
      considering pT sum of other particles within 
      DeltaR = 0.3. (No details given in the paper...) 

    Further requirements:
    - Identift the two highest pT EM objects (gg, ee, ge)
    - Require objects to have DeltaR > 0.6
    - Require objects to have m > 105

    Vetos:
    - Any *additional* electron with pT > 25, |eta| < 2.5
    - Any muon with pT > 25, |eta| < 2.4

    Six SRs based on MET value.    
*/

#include <vector>
#include <cmath>
#include <memory>
#include <iomanip>
#include <fstream>

#include "gambit/ColliderBit/analyses/Analysis.hpp"
#include "gambit/ColliderBit/CMSEfficiencies.hpp"
#include "gambit/ColliderBit/mt2_bisect.h"
#include "gambit/ColliderBit/analyses/Cutflow.hpp"

// #define CHECK_CUTFLOW

using namespace std;

namespace Gambit {
  namespace ColliderBit {


    class Analysis_CMS_13TeV_2Photon_GMSB_36invfb : public Analysis {
    public:

      static constexpr const char* detector = "CMS";

      // Counters for the number of accepted events for each signal region
      std::map<string,double> _numSR = {
        {"SR_MET_100-115",  0},
        {"SR_MET_115-130",  0},
        {"SR_MET_130-150",  0},
        {"SR_MET_150-185",  0},
        {"SR_MET_185-250",  0},
        {"SR_MET_>250",  0},
      };

      // Cutflow _cutflow;

      // Analysis_CMS_13TeV_2Photon_GMSB_36invfb():
      // _cutflow("CMS 2-photon GMSB 13 TeV", {"preselection", "MET>300GeV", "MT(g,MET)>300GeV", "S_T^g>600GeV"})
      Analysis_CMS_13TeV_2Photon_GMSB_36invfb()
      {
        set_analysis_name("CMS_13TeV_2Photon_GMSB_36invfb");
        set_luminosity(35.9);
      }


      void run(const HEPUtils::Event* event)
      {
        // Baseline objects
        HEPUtils::P4 ptot = event->missingmom();
        double met = event->met();

        // _cutflow.fillinit();

        // Photons
        // NOTE: 
        //   No photon efficiency info available for this analysis.
        //   We therefore assume the same efficiency map as used for 
        //   other CMS 36 fb^-1 SUSY searches:
        //
        //   https://twiki.cern.ch/twiki/pub/CMSPublic/SUSMoriond2017ObjectsEfficiency/PhotonEfficiencies_ForPublic_Moriond2017_LoosePixelVeto.pdf
        //
        //   The efficiency map has been extended to cover the low-pT region (pT < 20)
        const vector<double> aPhoton={0., 0.8, 1.4442, 1.566, 2.0, 2.5, DBL_MAX};  // Bin edges in eta
        const vector<double> bPhoton={0., 20., 35., 50., 90., DBL_MAX};  // Bin edges in pT. Assume flat efficiency above 500, where the CMS map stops.
        const vector<double> cPhoton={
                           // pT:   (0,20),  (20,35),  (35,50),  (50,90),  (90,inf)
                                     0.0,    0.735,    0.779,    0.805,    0.848,   // eta: (0, 0.8)
                                     0.0,    0.726,    0.746,    0.768,    0.809,   // eta: (0.8, 1.4442)
                                     0.0,    0.0,      0.0,      0.0,      0.0,     // eta: (1.4442, 1.566)
                                     0.0,    0.669,    0.687,    0.704,    0.723,   // eta: (1.566, 2.0)
                                     0.0,    0.564,    0.585,    0.592,    0.612,   // eta: (2.0, 2.5)
                                     0.0,    0.0,      0.0,      0.0,      0.0,     // eta > 2.5
                                 };
        HEPUtils::BinnedFn2D<double> _eff2dPhoton(aPhoton,bPhoton,cPhoton);
        vector<HEPUtils::Particle*> Photons;
        for (HEPUtils::Particle* photon : event->photons())
        {
          bool isPhoton=has_tag(_eff2dPhoton, photon->abseta(), photon->pT());
          if (isPhoton && photon->pT()>15.) Photons.push_back(photon);
        }


        // Electrons
        // NOTE:
        //   No electron efficiency info available for this analysis.
        //   We therefore assume the efficiency map used for the 36 fb^-1 CMS multilepton search:
        //   https://twiki.cern.ch/twiki/pub/CMSPublic/SUSMoriond2017ObjectsEfficiency/2d_full_pteta_el_039_multi_ttbar.pdf
        //
        //   See this page for more info:
        //   https://twiki.cern.ch/twiki/bin/view/CMSPublic/SUSMoriond2017ObjectsEfficiency
        //
        //   The efficiency map has been extended to cover the low-pT region, using the efficiencies from BuckFast (CMSEfficiencies.hpp)
        const vector<double> aEl={0., 0.8, 1.442, 1.556, 2., 2.5, DBL_MAX};   // Bin edges in eta
        const vector<double> bEl={0., 10., 15., 20., 25., 30., 40., 50., DBL_MAX}; // Bin edges in pT. Assume flat efficiency above 200, where the CMS map stops.
        const vector<double> cEl={
                          // pT: (0,10),  (10,15),  (15,20),  (20,25),  (25,30),  (30,40),  (40,50),  (50,inf)
                                   0.0,    0.95,    0.507,    0.619,    0.682,    0.742,    0.798,    0.863,  // eta: (0, 0.8)
                                   0.0,    0.95,    0.429,    0.546,    0.619,    0.710,    0.734,    0.833,  // eta: (0.8, 1.4429
                                   0.0,    0.95,    0.256,    0.221,    0.315,    0.351,    0.373,    0.437,  // eta: (1.442, 1.556)
                                   0.0,    0.85,    0.249,    0.404,    0.423,    0.561,    0.642,    0.749,  // eta: (1.556, 2)
                                   0.0,    0.85,    0.195,    0.245,    0.380,    0.441,    0.533,    0.644,  // eta: (2, 2.5)
                                   0.0,    0.0,     0.0,      0.0,      0.0,      0.0,      0.0,      0.0,    // eta > 2.5
                                  };
        HEPUtils::BinnedFn2D<double> _eff2dEl(aEl,bEl,cEl);
        vector<HEPUtils::Particle*> baselineElectrons;
        for (HEPUtils::Particle* electron : event->electrons()) {
          bool isEl=has_tag(_eff2dEl, electron->abseta(), electron->pT());
          // No info in the paper on pT or |eta| cuts for baseline electrons, 
          // but the above efficieny map effectively requires pT > 10 and |eta| < 2.5
          if (isEl) baselineElectrons.push_back(electron);
        }

        // Muons
        // NOTE:
        //   No muon efficiency info available for this analysis.
        //   We therefore assume the efficiency map used for the 36 fb^-1 CMS multilepton search:
        //   https://twiki.cern.ch/twiki/pub/CMSPublic/SUSMoriond2017ObjectsEfficiency/2d_full_pteta_mu_039_multi_ttbar.pdf
        //
        //   See this page for more info:
        //   https://twiki.cern.ch/twiki/bin/view/CMSPublic/SUSMoriond2017ObjectsEfficiency
        //
        //   The efficiency map has been extended to cover the low-pT region, using the efficiencies from BuckFast (CMSEfficiencies.hpp)
        const vector<double> aMu={0., 0.9, 1.2, 2.1, 2.4, DBL_MAX};   // Bin edges in eta
        const vector<double> bMu={0., 10., 15., 20., 25., 30., 40., 50., DBL_MAX};  // Bin edges in pT. Assume flat efficiency above 200, where the CMS map stops.
        const vector<double> cMu={
                           // pT:   (0,10),  (10,15),  (15,20),  (20,25),  (25,30),  (30,40),  (40,50),  (50,inf)
                                     0.0,     0.704,    0.797,    0.855,    0.880,    0.906,    0.927,    0.931,  // eta: (0, 0.9)
                                     0.0,     0.639,    0.776,    0.836,    0.875,    0.898,    0.940,    0.930,  // eta: (0.9, 1.2)
                                     0.0,     0.596,    0.715,    0.840,    0.862,    0.891,    0.906,    0.925,  // eta: (1.2, 2.1)
                                     0.0,     0.522,    0.720,    0.764,    0.803,    0.807,    0.885,    0.877,  // eta: (2.1, 2.4)
                                     0.0,     0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,    // eta > 2.4
                                 };
        HEPUtils::BinnedFn2D<double> _eff2dMu(aMu,bMu,cMu);
        vector<HEPUtils::Particle*> baselineMuons;
        for (HEPUtils::Particle* muon : event->muons()) {
          bool isMu=has_tag(_eff2dMu, muon->abseta(), muon->pT());
          // No info in the paper on pT or |eta| cuts for baseline muons, 
          // but the above efficieny map effectively requires pT > 10 and |eta| < 2.4
          if (isMu) baselineMuons.push_back(muon);
        }


        // _Anders: Got this far...


     //    // jets
     //    vector<HEPUtils::Jet*> Jets;
     //    for (HEPUtils::Jet* jet : event->jets())
     //    {
     //      if (jet->pT()>30. &&fabs(jet->eta())<3.0) Jets.push_back(jet);
     //    }
     //    // TODO: Apply jets isolation instead of removeOverlap.
     //    removeOverlap(Jets, Photons, 0.2);

     //    // Preselection
     //    bool high_pT_photon = false;  // At least one high-pT photon;
     //    bool delta_R_g_j = false;     // Photons are required to have delta_R>0.5 to the nearest jet;
     //    bool delta_phi_j_MET = false; // Jets with pT>100 GeV must fulfill delta_phi(MET,jet)>0.3;
      // for (HEPUtils::Particle* photon  : Photons){
      //     if (photon->pT()>180. && fabs(photon->eta()) < 1.44) {
      //         high_pT_photon = true;
      //         for (HEPUtils::Jet* jet : Jets){
      //             if ( jet->mom().deltaR_eta(photon->mom()) < 0.5 ) delta_R_g_j=true;
      //         }
      //     }
      // }
     //    if (not high_pT_photon) return;
     //    if (delta_R_g_j) return;
     //    for (HEPUtils::Jet* jet : Jets){
     //        if (jet->pT()>100. && jet->mom().deltaPhi(ptot) < 0.3 ) delta_phi_j_MET=true;
     //    }
     //    if (delta_phi_j_MET) return;

     //    // _cutflow.fill(1);


     //    // MET > 300 GeV
     //    if (met<300)return;
     //    // _cutflow.fill(2);

     //    // MT(photon,MET) > 300 GeV
     //    double MT = sqrt(2.*Photons[0]->pT()*met*(1. - std::cos(Photons[0]->mom().deltaPhi(ptot)) ));
     //    if (MT<300)return;
     //    // _cutflow.fill(3);

     //    // S_T^gamma > 600 GeV
     //    double STgamma = met;
     //    for (HEPUtils::Particle* photon  : Photons){
      //     STgamma += photon->pT();
      //   }
        
     //    if (STgamma<600) return;
     //    // _cutflow.fill(4);

     //    // Signal regions
     //    if      (STgamma<800)  _numSR["SR-600-800"]++;
     //    else if (STgamma<1000) _numSR["SR-800-1000"]++;
     //    else if (STgamma<1300) _numSR["SR-1000-1300"]++;
     //    else                   _numSR["SR-1300"]++;

      }


      /// Combine the variables of another copy of this analysis (typically on another thread) into this one.
      void combine(const Analysis* other)
      {
        const Analysis_CMS_13TeV_2Photon_GMSB_36invfb* specificOther
                = dynamic_cast<const Analysis_CMS_13TeV_2Photon_GMSB_36invfb*>(other);
        for (auto& el : _numSR) {
          el.second += specificOther->_numSR.at(el.first);
        }
      }


      virtual void collect_results()
      {
        // #ifdef CHECK_CUTFLOW
        // cout << _cutflow << endl;
        // for (auto& el : _numSR) {
        //     cout << el.first << "\t" << _numSR[el.first] << endl;
        // }
        // #endif

        // add_result(SignalRegionData("SR label", n_obs, {s, s_sys}, {b, b_sys}));
        add_result(SignalRegionData("SR_MET_100-115", 105, {_numSR["SR_MET_100-115"], 0.}, {114., 13.}));
        add_result(SignalRegionData("SR_MET_115-130", 39, {_numSR["SR_MET_115-130"], 0.}, {42.9, 7.5}));
        add_result(SignalRegionData("SR_MET_130-150", 21, {_numSR["SR_MET_130-150"], 0.}, {27.3, 5.6}));
        add_result(SignalRegionData("SR_MET_150-185", 21, {_numSR["SR_MET_150-185"], 0.}, {17.4, 4.1}));
        add_result(SignalRegionData("SR_MET_185-250", 11, {_numSR["SR_MET_185-250"], 0.}, {10.2, 2.7}));
        add_result(SignalRegionData("SR_MET_>250", 12, {_numSR["SR_MET_>250"], 0.}, {5.4, 1.6}));
      }


    protected:
      void analysis_specific_reset() {
       for (auto& el : _numSR) { el.second = 0.;}
      }

    };

    // Factory fn
    DEFINE_ANALYSIS_FACTORY(CMS_13TeV_2Photon_GMSB_36invfb)


  }
}
