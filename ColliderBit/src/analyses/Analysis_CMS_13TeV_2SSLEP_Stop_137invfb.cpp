///
///  \author Yang Zhang
///  \date 2019 June
///  *********************************************

// Based on http://cms-results.web.cern.ch/cms-results/public-results/preliminary-results/SUS-19-008/index.html
// Search for physics beyond the standard model in events with two same-sign leptons or at least three leptons and jets in proton-proton collisions at 13 TeV

// Note:
// 1. Not fully validated.

#include <vector>
#include <cmath>
#include <memory>
#include <iomanip>

#include "gambit/ColliderBit/analyses/Analysis.hpp"
#include "gambit/ColliderBit/mt2_bisect.h"
#include "gambit/ColliderBit/CMSEfficiencies.hpp"
#include "gambit/ColliderBit/analyses/Cutflow.hpp"

// #define CHECK_CUTFLOW

using namespace std;

namespace Gambit {
  namespace ColliderBit {

    class Analysis_CMS_13TeV_2SSLEP_Stop_137invfb : public Analysis {
    protected:
        // Counters for the number of accepted events for each signal region
        static const size_t NUMSRHH = 62;
        double _SRHH[NUMSRHH];
        static const size_t NUMSRHL = 43;
        double _SRHL[NUMSRHL];
        static const size_t NUMSRLL = 8;
        double _SRLL[NUMSRLL];
        static const size_t NUMSRLM = 11;
        double _SRLM[NUMSRLM];
        static const size_t NUMSRML = 44;
        double _SRML[NUMSRML];
        Cutflow _cutflow;


    public:

        // Required detector sim
        static constexpr const char* detector = "CMS";

        Analysis_CMS_13TeV_2SSLEP_Stop_137invfb():
        _cutflow("CMS_13TeV_2SSLEP_Stop_137invfb", {"Trigger_and_2leptons", "At_least_one_SS_lepton_pair", "Baseline"})
        {

            set_analysis_name("CMS_13TeV_2SSLEP_Stop_137invfb");
            set_luminosity(137);
            
            for (size_t i = 0; i < NUMSRHH; ++i) _SRHH[i] = 0;
            for (size_t i = 0; i < NUMSRHL; ++i) _SRHL[i] = 0;
            for (size_t i = 0; i < NUMSRLL; ++i) _SRLL[i] = 0;
            for (size_t i = 0; i < NUMSRLM; ++i) _SRLM[i] = 0;
            for (size_t i = 0; i < NUMSRML; ++i) _SRML[i] = 0;
        }

        struct ptComparison {
            bool operator() (HEPUtils::Particle* i,HEPUtils::Particle* j) {return (i->pT()>j->pT());}
        } comparePt;
        
        void run(const HEPUtils::Event* event) {
            _cutflow.fillinit();
            
            // Missing energy
            double met = event->met();
            HEPUtils::P4 ptot = event->missingmom();

            // Electrons
            //@note Numbers digitized from https://twiki.cern.ch/twiki/pub/CMSPublic/SUSMoriond2017ObjectsEfficiency/2d_full_pteta_el_035_ttbar.pdf
            const vector<double> aEl={0., 0.8, 1.442, 1.556, 2., 2.5, DBL_MAX};   // Bin edges in eta
            const vector<double> bEl={0., 15., 20., 25., 30., 40., 50, DBL_MAX}; // Bin edges in pT. Assume flat efficiency above 200, where the CMS map stops.
            const vector<double> cEl={
                          // pT:  (0,15), (15,20), (20,25), (25,30), (30,40), (40,50), (50,inf)
                                   0.0,   0.398,   0.501,   0.556,   0.619,   0.669,   0.720,// eta: (0, 0.8)
                                   0.0,   0.344,   0.433,   0.498,   0.579,   0.600,   0.671,// eta: (0.8, 1.4429)
                                   0.0,   0.201,   0.156,   0.206,   0.222,   0.255,   0.307,// eta: (1.442, 1.556)
                                   0.0,   0.210,   0.302,   0.338,   0.428,   0.484,   0.561,// eta: (1.556, 2)
                                   0.0,   0.162,   0.172,   0.250,   0.339,   0.396,   0.444,// eta: (2, 2.5)
                                   0.0,   0.0,     0.0,     0.0,     0.0,     0.0,     0.0// eta > 2.5
                                  };
            HEPUtils::BinnedFn2D<double> _eff2dEl(aEl,bEl,cEl);
            vector<HEPUtils::Particle*> electrons;
            for (HEPUtils::Particle* electron : event->electrons()) {
                bool isEl=has_tag(_eff2dEl, fabs(electron->eta()), electron->pT());
                if (electron->pT() > 15. && fabs(electron->eta()) < 2.5 && isEl)
                    electrons.push_back(electron);
            }
            
            // Muons
            //@note Numbers digitized from https://twiki.cern.ch/twiki/pub/CMSPublic/SUSMoriond2017ObjectsEfficiency/2d_full_pteta_mu_035_ttbar.pdf
            const vector<double> aMu={0., 0.9, 1.2, 2.1, 2.4, DBL_MAX};   // Bin edges in eta
            const vector<double> bMu={0., 10, 15., 20., 25., 30, 40, 50, DBL_MAX};  // Bin edges in pT. Assume flat efficiency above 200, where the CMS map stops.
            const vector<double> cMu={
                          // pT:  (0,10), (10,15), (15,20), (20,25), (25,30), (30,40), (40,50), (50,inf)
                                   0.0,   0.564,   0.645,    0.739,  0.803,   0.860,   0.894,   0.907, // eta: (0, 0.9)
                                   0.0,   0.525,   0.616,    0.700,  0.773,   0.825,   0.891,   0.898, // eta: (0.9, 1.2)
                                   0.0,   0.514,   0.572,    0.697,  0.748,   0.789,   0.837,   0.870, // eta: (1.2, 2.1)
                                   0.0,   0.440,   0.575,    0.604,  0.663,   0.696,   0.784,   0.794,// eta: (2.1, 2.4)
                                   0.0,   0.0,     0.0,      0.0,    0.0,     0.0,     0.0,     0.0// eta > 2.4
                                  };
            HEPUtils::BinnedFn2D<double> _eff2dMu(aMu,bMu,cMu);
            vector<HEPUtils::Particle*> muons;
            for (HEPUtils::Particle* muon : event->muons()) {
                bool isMu=has_tag(_eff2dMu, fabs(muon->eta()), muon->pT());
                if (muon->pT() > 10.&& fabs(muon->eta()) < 2.4 && isMu)
                    muons.push_back(muon);
            }
            
            double HT;
            // Jets
            vector<HEPUtils::Jet*> candJets;
            for (HEPUtils::Jet* jet : event->jets()) {
                if (jet->pT() > 25. && fabs(jet->eta()) < 2.4){
                    HT += jet->pT();
                    candJets.push_back(jet);
                }
            }

            // Jets
            vector<HEPUtils::Jet*> bJets;
            vector<HEPUtils::Jet*> nonbJets;
            
            // Find b-jets
            // Copied from ATLAS_13TeV_3b_24invfb
            double btag = 0.85; double cmisstag = 1/12.; double misstag = 1./381.;
            for (HEPUtils::Jet* jet : candJets) {
                // Tag
                if( jet->btag() && random_bool(btag) ) bJets.push_back(jet);
                // Misstag c-jet
                else if( jet->ctag() && random_bool(cmisstag) ) bJets.push_back(jet);
                // Misstag light jet
                else if( random_bool(misstag) ) bJets.push_back(jet);
                // Non b-jet
                else if( jet->pT() > 40. ) nonbJets.push_back(jet);
            }

//            // Overlap removal
//            JetLeptonOverlapRemoval(candJets,electrons,0.2);
//            LeptonJetOverlapRemoval(electrons,candJets);
//            JetLeptonOverlapRemoval(candJets,muons,0.4);
//            LeptonJetOverlapRemoval(muons,candJets);
            
            size_t Nb=bJets.size();
            size_t Nj=nonbJets.size();
            
            // Leptons = electrons + muons
            vector<HEPUtils::Particle*> leptons;
            leptons=electrons;
            leptons.insert(leptons.end(),muons.begin(),muons.end());
            sort(leptons.begin(),leptons.end(),comparePt);

            // At least two light leptons
            if (leptons.size()<2) return;
            
            // Find pair same sign (SS) leptons
            vector<size_t> SS_1,SS_2;
            for (size_t i=0; i<leptons.size(); ++i) {
                for (size_t j=i+1; j<leptons.size(); ++j) {
                    if (leptons[i]->pid()*leptons[j]->pid()>0){
                        SS_1.push_back(i);
                        SS_2.push_back(j);
                    }
                    // mll>12 for an opposite-sign same flavor pair lepton
                    if (leptons[i]->pid()+leptons[j]->pid()==0 and (leptons[i]->mom()+leptons[j]->mom()).m()<12) return;
                    // mll>8 GeV for any pair of leptons
                    if ((leptons[i]->mom()+leptons[j]->mom()).m()<8) return;
                }
            }
            _cutflow.fill(1);
            
            // One SS lepton pair
            if (SS_1.size()==0) return;
            _cutflow.fill(2);
            
            // At least two jets and MET>50
            if (nonbJets.size()<2 or  met<50) return;
            _cutflow.fill(3);
            
            // Find the only SS lepton pair
            size_t SS1 = SS_1[0];
            size_t SS2 = SS_2[0];
            bool find_one_muon = false;
            for (size_t i=1; i<SS_1.size(); ++i) {
                // SS_1 and SS_2 are already order by lepton PT sum
                if (fabs(leptons[SS_1[i]]->pid())==13 and fabs(leptons[SS_1[i]]->pid())==13) {
                    // both of the leptons are muon
                    SS1 = SS_1[i];
                    SS2 = SS_2[i];
                    break;
                }
                if ( (not find_one_muon) and (fabs(leptons[SS_1[i]]->pid())==13 or fabs(leptons[SS_1[i]]->pid())==13)){
                    // one of the leptons is muon
                    SS1 = SS_1[i];
                    SS2 = SS_2[i];
                    find_one_muon = true;
                }
            }
            
            // M_T^{miss}
            double MTmiss = sqrt(2.*leptons[SS1]->pT()*met*(1. - cos(leptons[SS1]->mom().deltaPhi(ptot))));
            if (MTmiss>sqrt(2.*leptons[SS2]->pT()*met*(1. - cos(leptons[SS2]->mom().deltaPhi(ptot))))) {
                MTmiss = sqrt(2.*leptons[SS2]->pT()*met*(1. - cos(leptons[SS2]->mom().deltaPhi(ptot))));
            }
            
            bool pp = leptons[SS1]->pid()>0;
            bool met_50_200 = met>50 and met<200;
            bool met_200_300 = met>200 and met<300;
            bool met_300_500 = met>300 and met<500;
            bool met_500 = met>500;
            bool MTmiss_l_120 = MTmiss<120;
            bool MTmiss_g_120 = MTmiss>120;
            bool Nj_2_4 = Nj>=2 and Nj<=4;
            bool Nj_5 = Nj>=5;
            bool HT_300 = HT<300;
            bool HT_300_1125 = HT>300 and HT<1125;
            bool HT_1125_1300 = HT>1125 and HT<1300;
            bool HT_1300_1600 = HT>1300 and HT<1600;
            bool HT_1300 = HT>1300;
            bool HT_1600 = HT>1600;
            bool SSHH_combine = ((MTmiss_l_120&&((met_50_200&&Nj_5)||met_200_300)) or (MTmiss_g_120&&met<300) )&&HT_300 ;
            
            // SSHH: exactly 2 leptons, both with PT>25 GeV, and MET>50 GeV
            if (leptons.size()==2 and leptons[1]->pT() > 25.) {
                if (Nb==0) {
                    if (MTmiss_l_120 and met_50_200  and Nj_2_4 and HT_300)               _SRHH[0]++;
                    if (MTmiss_l_120 and met_50_200  and Nj_2_4 and HT_300_1125)          _SRHH[1]++;
                    if (SSHH_combine)                                                     _SRHH[2]++;
                    if (MTmiss_l_120 and met_50_200  and Nj_5   and HT_300_1125)          _SRHH[3]++;
                    if (MTmiss_l_120 and met_200_300 and Nj_2_4 and HT_300_1125 and pp)   _SRHH[4]++;
                    if (MTmiss_l_120 and met_200_300 and Nj_2_4 and HT_300_1125 and !pp)  _SRHH[5]++;
                    if (MTmiss_l_120 and met_200_300 and Nj_5   and HT_300_1125)          _SRHH[6]++;
                    if (MTmiss_g_120 and met_50_200  and Nj_2_4 and HT_300_1125 and pp)   _SRHH[7]++;
                    if (MTmiss_g_120 and met_50_200  and Nj_2_4 and HT_300_1125 and !pp)  _SRHH[8]++;
                    if (MTmiss_g_120 and((met_50_200&&Nj_5)||met_200_300)and HT_300_1125) _SRHH[9]++;
                } else if (Nb==1) {
                    if (MTmiss_l_120 and met_50_200  and Nj_2_4 and HT_300)              _SRHH[10]++;
                    if (MTmiss_l_120 and met_50_200  and Nj_2_4 and HT_300_1125)         _SRHH[11]++;
                    if (SSHH_combine and pp)                                             _SRHH[12]++;
                    if (SSHH_combine and !pp)                                            _SRHH[13]++;
                    if (MTmiss_l_120 and met_50_200  and Nj_5   and HT_300_1125 and pp)  _SRHH[14]++;
                    if (MTmiss_l_120 and met_50_200  and Nj_5   and HT_300_1125 and !pp) _SRHH[15]++;
                    if (MTmiss_l_120 and met_200_300 and Nj_2_4 and HT_300_1125 and pp)  _SRHH[16]++;
                    if (MTmiss_l_120 and met_200_300 and Nj_2_4 and HT_300_1125 and !pp) _SRHH[17]++;
                    if (MTmiss_l_120 and met_200_300 and Nj_5   and HT_300_1125)         _SRHH[18]++;
                    if (MTmiss_g_120 and met_50_200  and Nj_2_4 and HT_300_1125 and pp)  _SRHH[19]++;
                    if (MTmiss_g_120 and met_50_200  and Nj_2_4 and HT_300_1125 and !pp) _SRHH[20]++;
                    if (MTmiss_g_120 and((met_50_200&&Nj_5)||met_200_300)and HT_300_1125) _SRHH[21]++;
                } else if (Nb==2){
                    if (MTmiss_l_120 and met_50_200  and Nj_2_4 and HT_300)              _SRHH[22]++;
                    if (MTmiss_l_120 and met_50_200  and Nj_2_4 and HT_300_1125)         _SRHH[23]++;
                    if (SSHH_combine and pp)                                             _SRHH[24]++;
                    if (SSHH_combine and !pp)                                            _SRHH[25]++;
                    if (MTmiss_l_120 and met_50_200  and Nj_5   and HT_300_1125 and pp)  _SRHH[26]++;
                    if (MTmiss_l_120 and met_50_200  and Nj_5   and HT_300_1125 and !pp) _SRHH[27]++;
                    if (MTmiss_l_120 and met_200_300 and Nj_2_4 and HT_300_1125 and pp)  _SRHH[28]++;
                    if (MTmiss_l_120 and met_200_300 and Nj_2_4 and HT_300_1125 and !pp) _SRHH[29]++;
                    if (MTmiss_l_120 and met_200_300 and Nj_5   and HT_300_1125)         _SRHH[30]++;
                    if (MTmiss_g_120 and met_50_200  and Nj_2_4 and HT_300_1125 and pp)  _SRHH[31]++;
                    if (MTmiss_g_120 and met_50_200  and Nj_2_4 and HT_300_1125 and !pp) _SRHH[32]++;
                    if (MTmiss_g_120 and((met_50_200&&Nj_5)||met_200_300)and HT_300_1125) _SRHH[33]++;
                } else if (Nb>=3){
                    if (MTmiss_l_120 and met<300                and HT_300 and pp)       _SRHH[34]++;
                    if (MTmiss_l_120 and met<300                and HT_300 and !pp)      _SRHH[35]++;
                    if (MTmiss_l_120 and met<300     and Nj_2_4 and HT_300_1125 and pp)  _SRHH[36]++;
                    if (MTmiss_l_120 and met<300     and Nj_2_4 and HT_300_1125 and !pp) _SRHH[37]++;
                    if (MTmiss_l_120 and met<300     and Nj_5   and HT_300_1125 and pp)  _SRHH[38]++;
                    if (MTmiss_l_120 and met<300     and Nj_5   and HT_300_1125 and !pp) _SRHH[39]++;
                    if (MTmiss_g_120                            and HT_300)              _SRHH[40]++;
                    if (MTmiss_g_120 and met<300     and Nj_2_4 and HT_300_1125 and pp)  _SRHH[41]++;
                    if (MTmiss_g_120 and met<300     and Nj_2_4 and HT_300_1125 and !pp) _SRHH[42]++;
                    if (MTmiss_g_120 and met<300     and Nj_5   and HT_300_1125 and pp)  _SRHH[43]++;
                    if (MTmiss_g_120 and met<300     and Nj_5   and HT_300_1125 and !pp) _SRHH[44]++;
                }
                
                if (met_300_500 and Nj_2_4 and HT>300 and pp)  _SRHH[45]++;
                if (met_300_500 and Nj_2_4 and HT>300 and !pp) _SRHH[46]++;
                if (met_500     and Nj_2_4 and HT>300 and pp)  _SRHH[47]++;
                if (met_500     and Nj_2_4 and HT>300 and !pp) _SRHH[48]++;
                if (met_300_500 and Nj_5   and HT>300 and pp)  _SRHH[49]++;
                if (met_300_500 and Nj_5   and HT>300 and !pp) _SRHH[50]++;
                if (met_500     and Nj_5   and HT>300 and pp)  _SRHH[51]++;
                if (met_500     and Nj_5   and HT>300 and !pp) _SRHH[52]++;
                
                if (HT_1125_1300 and met<300 and Nj<5) _SRHH[53]++;
                if (HT_1300_1600 and met<300 and Nj<5) _SRHH[54]++;
                if (HT_1600      and met<300 and Nj<5) _SRHH[55]++;
                
                if (HT_1125_1300 and met<300 and Nj==5 and Nj==6) _SRHH[56]++;
                if (HT_1300_1600 and met<300 and Nj==5 and Nj==6) _SRHH[57]++;
                if (HT_1600      and met<300 and Nj==5 and Nj==6) _SRHH[58]++;
                
                if (HT_1125_1300 and met<300 and Nj<5) _SRHH[59]++;
                if (HT_1300_1600 and met<300 and Nj<5) _SRHH[60]++;
                if (HT_1600      and met<300 and Nj<5) _SRHH[61]++;
            }
            
            bool SSHL_combine = MTmiss_l_120&&( (met_50_200&&Nj_5) or met_200_300 )&&HT_300 ;
            
            // SSHL: exactly 2 leptons, one with PT>25 GeV, one with PT<25 GeV,  and MET>50 GeV
            if (leptons.size()==2 and leptons[0]->pT() > 25. and leptons[1]->pT() < 25.) {
                if (Nb==0 and MTmiss_l_120) {
                    if ( met_50_200  and Nj_2_4 and HT_300)              _SRHL[0]++;
                    if ( met_50_200  and Nj_2_4 and HT_300_1125)         _SRHL[1]++;
                    if ( SSHL_combine)                                   _SRHL[2]++;
                    if ( met_50_200  and Nj_5 and HT_300_1125)           _SRHL[3]++;
                    if ( met_200_300 and Nj_2_4 and HT_300_1125 and pp)  _SRHL[4]++;
                    if ( met_200_300 and Nj_2_4 and HT_300_1125 and !pp) _SRHL[5]++;
                    if ( met_200_300 and Nj_5 and HT_300_1125)           _SRHL[6]++;
                    
                } else if(Nb==1 and MTmiss_l_120) {
                    if ( met_50_200  and Nj_2_4 and HT_300)              _SRHL[7]++;
                    if ( met_50_200  and Nj_2_4 and HT_300_1125)         _SRHL[8]++;
                    if ( SSHL_combine and pp)                            _SRHL[9]++;
                    if ( SSHL_combine and !pp)                           _SRHL[10]++;
                    if ( met_50_200  and Nj_5 and HT_300_1125 and pp)    _SRHL[11]++;
                    if ( met_50_200  and Nj_5 and HT_300_1125 and !pp)   _SRHL[12]++;
                    if ( met_200_300 and Nj_2_4 and HT_300_1125)         _SRHL[13]++;
                    if ( met_200_300 and Nj_5 and HT_300_1125 and pp)    _SRHL[14]++;
                    if ( met_200_300 and Nj_5 and HT_300_1125 and !pp)   _SRHL[15]++;
                } else if(Nb==2 and MTmiss_l_120) {
                    if ( met_50_200  and Nj_2_4 and HT_300)              _SRHL[16]++;
                    if ( met_50_200  and Nj_2_4 and HT_300_1125)         _SRHL[17]++;
                    if ( SSHL_combine and pp)                            _SRHL[18]++;
                    if ( SSHL_combine and !pp)                           _SRHL[19]++;
                    if ( met_50_200  and Nj_5 and HT_300_1125 and pp)    _SRHL[20]++;
                    if ( met_50_200  and Nj_5 and HT_300_1125 and !pp)   _SRHL[21]++;
                    if ( met_200_300 and Nj_2_4 and HT_300_1125 and pp)  _SRHL[22]++;
                    if ( met_200_300 and Nj_2_4 and HT_300_1125 and !pp) _SRHL[23]++;
                    if ( met_200_300 and Nj_5 and HT_300_1125)           _SRHL[24]++;
                }else if(Nb==3 and MTmiss_l_120) {
                    if ( met_50_200 and HT_300 and pp)                   _SRHL[25]++;
                    if ( met_50_200 and HT_300 and !pp)                  _SRHL[26]++;
                    if ( met_50_200 and HT_300_1125 and pp)              _SRHL[27]++;
                    if ( met_50_200 and HT_300_1125 and !pp)             _SRHL[28]++;
                    if ( met_200_300 and HT_300_1125)                    _SRHL[29]++;
                }
                if (MTmiss_g_120 and met<300 and HT_300)       _SRHL[30]++;
                if (MTmiss_g_120 and met<300 and HT_300_1125)  _SRHL[31]++;
                if (met_300_500 and Nj_2_4 and HT>300 and pp)  _SRHL[32]++;
                if (met_300_500 and Nj_2_4 and HT>300 and !pp) _SRHL[33]++;
                if (met_500 and Nj_2_4 and HT>300 and pp)      _SRHL[34]++;
                if (met_500 and Nj_2_4 and HT>300 and !pp)     _SRHL[35]++;
                if (met_300_500 and Nj_5 and HT>300 and pp)    _SRHL[36]++;
                if (met_300_500 and Nj_5 and HT>300 and !pp)   _SRHL[37]++;
                if (met_500 and Nj_5 and HT>300)               _SRHL[38]++;
                
                if (HT_1125_1300 and pp)  _SRHL[39]++;
                if (HT_1125_1300 and !pp) _SRHL[40]++;
                if (HT_1300 and pp)       _SRHL[41]++;
                if (HT_1300 and !pp)      _SRHL[42]++;
            }

            // SSLL: exactly 2 leptons, both with PT<25 GeV, and MET>50 GeV
            if (leptons.size()==2 and leptons[0]->pT() < 25. and leptons[1]->pT() < 25.) {
                if (HT>400) {
                    if (MTmiss_l_120) {
                        if (Nb==0) {
                            if (met_50_200) _SRLL[0]++;
                            else            _SRLL[1]++;
                        } else if (Nb==1) {
                            if (met_50_200) _SRLL[2]++;
                            else            _SRLL[3]++;
                        } else if (Nb==2) {
                            if (met_50_200) _SRLL[4]++;
                            else            _SRLL[5]++;
                        } else if (Nb>=3)   _SRLL[6]++;
                    } else                  _SRLL[7]++;
                }
            }

            // LM: exactly 2 leptons, both with PT>25 GeV, and MET<50 GeV
            if (leptons.size()==2 and leptons[0]->pT() > 25. and leptons[1]->pT() > 25.) {
                if (HT_300_1125) {
                    if (Nb==0) {
                        if (Nj_2_4) _SRLM[0]++;
                        if (Nj_5)   _SRLM[1]++;
                    } else if (Nb==1){
                        if (Nj_2_4) _SRLM[2]++;
                        if (Nj_5)   _SRLM[3]++;
                    } else if (Nb==2){
                        if (Nj_2_4) _SRLM[4]++;
                        if (Nj_5)   _SRLM[5]++;
                    } else          _SRLM[6]++;
                } else if (HT_1125_1300){
                    if (Nj_2_4) {
                        _SRLM[7]++;
                    } else if (Nj_5) {
                        _SRLM[8]++;
                    }
                } else if (HT_1300){
                    if (Nj_2_4) {
                        _SRLM[9]++;
                    } else if (Nj_5) {
                        _SRLM[10]++;
                    }
                }
            }
            
            // ML: >=3 leptons, at least one with PT>25 GeV, and MET>50 GeV
            if (leptons.size()>=3 and leptons[0]->pT() > 25.) {
                bool on_Z=false;
                size_t Z_l1;
                size_t Z_l2;
                for (size_t i=0; i<leptons.size(); ++i) {
                    for (size_t j=i+1; j<leptons.size(); ++j) {
                        if (leptons[i]->pid()+leptons[j]->pid()==0 and fabs((leptons[i]->mom()+leptons[j]->mom()).m()-91)<15) {
                            on_Z = true;
                            Z_l1=i;
                            Z_l2=j;
                            break;
                        }
                    }
                    if (on_Z) break;
                }
                                
                if (not on_Z){ // off-Z
                    if (Nb==0) {
                        if (HT<400) {
                            if (met<150) {
                                if (MTmiss_l_120) _SRML[0]++;
                                else              _SRML[1]++;
                            } else if (met<300) {
                                if (MTmiss_l_120) _SRML[2]++;
                                else              _SRML[3]++;
                            }
                        }
                        else if (HT<600) {
                             if (met<150)         _SRML[4]++;
                             else if (met<300)    _SRML[5]++;
                        }
                    } else if (Nb==1) {
                        if (HT<400) {
                             if (met<150)         _SRML[6]++;
                             else if (met<300)    _SRML[7]++;
                        }
                        else if (HT<600) {
                             if (met<150)         _SRML[8]++;
                             else if (met<300)    _SRML[9]++;
                        }
                    } else if (Nb==2) {
                        if (HT<400) {
                             if (met<150)         _SRML[10]++;
                             else if (met<300)    _SRML[11]++;
                        }
                        else if (HT<600) {
                             if (met<150)         _SRML[12]++;
                             else if (met<300)    _SRML[13]++;
                        }
                    } else if (Nb>=3) {
                        if (HT<600)               _SRML[14]++;
                    }
                    if (HT>=600) {
                        if (met<150) {
                            if (MTmiss_l_120)     _SRML[15]++;
                            else                  _SRML[16]++;
                        } else if (met<300) {
                            if (MTmiss_l_120)     _SRML[17]++;
                            else                  _SRML[18]++;
                        }
                    }
                    if (met>300) {
                        if (MTmiss_l_120)         _SRML[19]++;
                        else                      _SRML[20]++;
                    }
                } else { // on-Z
                    // M_T^{miss} for on-Z
                    double MTmiss_Z = 1000;  
                    for (size_t i=0; i<leptons.size(); ++i) {
                        if ( (i != Z_l1) and (i != Z_l2)) {
                            double MTmiss_try = sqrt(2.*leptons[i]->pT()*met*(1. - cos(leptons[i]->mom().deltaPhi(ptot))));
                            if (MTmiss_try<MTmiss_Z) MTmiss_Z=MTmiss_try;
                        }
                    }
                    bool MTmiss_Z_l_120 = MTmiss_Z<120;
                    if (Nb==0) {
                        if (HT<400) {
                            if (met<150) {
                                if (MTmiss_Z_l_120) _SRML[21]++;
                                else                _SRML[22]++;
                            } else if (met<300) {
                                if (MTmiss_Z_l_120) _SRML[23]++;
                                else                _SRML[24]++;
                            }
                        }
                        else if (HT<600) {
                            if (met<150) {
                                if (MTmiss_Z_l_120) _SRML[25]++;
                                else                _SRML[26]++;
                            } else if (met<300) {
                                if (MTmiss_Z_l_120) _SRML[27]++;
                                else                _SRML[28]++;
                            }
                        }
                    } else if (Nb==1) {
                        if (HT<400) {
                             if (met<150)           _SRML[29]++;
                             else if (met<300)      _SRML[30]++;
                        }
                        else if (HT<600) {
                             if (met<150)           _SRML[31]++;
                             else if (met<300)      _SRML[32]++;
                        }
                    } else if (Nb==2) {
                        if (HT<400) {
                             if (met<150)           _SRML[33]++;
                             else if (met<300)      _SRML[34]++;
                        }
                        else if (HT<600) {
                             if (met<150)           _SRML[35]++;
                             else if (met<300)      _SRML[36]++;
                        }
                    } else if (Nb>=3) {
                        if (HT<600)                 _SRML[37]++;
                    }
                    if (HT>=600) {
                        if (met<150) {
                            if (MTmiss_Z_l_120)     _SRML[38]++;
                            else                    _SRML[39]++;
                        } else if (met<300) {
                            if (MTmiss_Z_l_120)     _SRML[40]++;
                            else                    _SRML[41]++;
                        }
                    }
                    if (met>300) {
                        if (MTmiss_Z_l_120)         _SRML[42]++;
                        else                        _SRML[43]++;
                    }
                }
            }
            
            return;
        }

        /// Combine the variables of another copy of this analysis (typically on another thread) into this one.
        void combine(const Analysis* other)
        {
            const Analysis_CMS_13TeV_2SSLEP_Stop_137invfb* specificOther
                = dynamic_cast<const Analysis_CMS_13TeV_2SSLEP_Stop_137invfb*>(other);
            
            for (size_t i = 0; i < NUMSRHH; ++i) _SRHH[i] += specificOther->_SRHH[i];
            for (size_t i = 0; i < NUMSRHL; ++i) _SRHL[i] += specificOther->_SRHL[i];
            for (size_t i = 0; i < NUMSRLL; ++i) _SRLL[i] += specificOther->_SRLL[i];
            for (size_t i = 0; i < NUMSRLM; ++i) _SRLM[i] += specificOther->_SRLM[i];
            for (size_t i = 0; i < NUMSRML; ++i) _SRML[i] += specificOther->_SRML[i];
        }


        void collect_results() {

            // Observed event counts
            static const double OBSNUM_HH[NUMSRHH] = {
                1510, 590, 103, 38, 57, 32, 5.5, 25, 21, 9.4, 930, 330, 36, 25, 44, 39, 27, 14.8, 11.5, 11.8, 9.6, 10.0, 270, 143, 15.2, 13.8, 33, 29, 11.5, 6.7, 7.5, 5.9, 6.5, 6.7, 10.3, 8.6, 10.6, 7.3, 9.6, 9.2, 1.3, 0.6, 0.8, 0.7, 0.7, 42, 18, 13, 2.0, 6.3, 3.7, 1.26, 0.4, 10.1, 7.0, 4.3, 5.3, 6, 2.2, 1.8, 1.9, 1.3
            };
            // Background estimates
            static const double BKGNUM_HH[NUMSRHH] = {
                310, 90, 22, 7, 10, 9, 1.7, 6, 5, 1.9, 230, 70, 7, 5, 7, 8, 5, 3.2, 3.0, 2.6, 2.1, 1.6, 40, 20, 2.4, 3.4, 5, 4, 2.5, 1.2, 1.8, 1.0, 1.9, 1.2, 1.9, 1.7, 2.0, 1.3, 2.2, 1.9, 0.6, 0.4, 0.4, 0.4, 0.5, 7, 4, 9, 0.5, 1.0, 0.7, 0.33, 0.4, 1.5, 1.1, 0.9, 0.8, 6, 0.4, 0.5, 0.4, 0.9
            };
            // Background uncertainties
            static const double BKGERR_HH[NUMSRHH] = {
                1609, 647, 132, 51, 49, 23, 7, 31, 20, 11, 1068, 370, 38, 31, 63, 38, 30, 15, 12, 14, 16, 15, 345, 169, 11, 18, 43, 38, 9, 5, 6, 14, 7, 11, 17, 11, 6, 5, 8, 11, 2, 1, 0, 1, 1, 59, 23, 10, 4, 13, 4, 4, 2, 24, 4, 5, 7, 6, 3, 5, 4, 0
            };
            for (size_t ibin = 0; ibin < NUMSRHH; ++ibin)
            {
                stringstream ss; ss << "SR-HH-" << ibin+1;
                add_result(SignalRegionData(ss.str(), OBSNUM_HH[ibin], {_SRHH[ibin],  0.}, {BKGNUM_HH[ibin], BKGERR_HH[ibin]}));
            }

            // Observed event counts
            static const double OBSNUM_HL[NUMSRHL] = {
                1300, 310, 25, 32, 29, 17, 4.5, 1010, 270, 7.1, 6.5, 39, 31, 23, 2.1, 1.7, 210, 85, 2.5, 3.0, 18.9, 15.9, 3.3, 4.4, 4.5, 8.2, 8.1, 9.7, 10.8, 1.1, 2.2, 2.6, 22, 7.2, 2.3, 0.42, 3.2, 1.4, 0.41, 3.1, 4, 4.7, 1.71
            };
            // Background estimates
            static const double BKGNUM_HL[NUMSRHL] = {
                310, 70, 6, 8, 6, 6, 2.5, 250, 60, 1.7, 1.6, 9, 8, 5, 1.1, 0.9, 40, 14, 1.2, 1.5, 3.5, 2.8, 0.6, 1.6, 1.7, 2.2, 2.2, 2.1, 2.8, 0.4, 0.5, 0.5, 6, 1.4, 0.5, 0.33, 1.5, 0.6, 0.25, 0.7, 4, 0.9, 0.35
            };
            // Background uncertainties
            static const double BKGERR_HL[NUMSRHL] = {
                1504, 319, 32, 32, 32, 11, 6, 1223, 307, 5, 7, 42, 37, 27, 7, 2, 256, 104, 4, 3, 27, 18, 2, 2, 5, 8, 6, 12, 7, 3, 5, 3, 23, 8, 4, 1, 3, 0, 0, 7, 0, 8, 6
            };
            for (size_t ibin = 0; ibin < NUMSRHL; ++ibin)
            {
                stringstream ss; ss << "SR-HL-" << ibin+1;
                add_result(SignalRegionData(ss.str(), OBSNUM_HL[ibin], {_SRHL[ibin],  0.}, {BKGNUM_HL[ibin], BKGERR_HL[ibin]}));
            }

            // Observed event counts
            static const double OBSNUM_LL[NUMSRLL] = {
                20, 4.9, 23, 4, 7.7, 1.5, 2.5, 0.07
            };
            // Background estimates
            static const double BKGNUM_LL[NUMSRLL] = {
                6, 1.5, 6, 4, 2.1, 0.8, 1.2, 0.07
            };
            // Background uncertainties
            static const double BKGERR_LL[NUMSRLL] = {
                25, 6, 29, 8, 13, 0, 0, 0
            };
            for (size_t ibin = 0; ibin < NUMSRLL; ++ibin)
            {
                stringstream ss; ss << "SR-LL-" << ibin+1;
                add_result(SignalRegionData(ss.str(), OBSNUM_LL[ibin], {_SRLL[ibin],  0.}, {BKGNUM_LL[ibin], BKGERR_LL[ibin]}));
            }

            // Observed event counts
            static const double OBSNUM_LM[NUMSRLM] = {
                240, 20, 140, 32, 54, 22, 9.7, 1.5, 1.6, 2.9, 1.9
            };
            // Background estimates
            static const double BKGNUM_LM[NUMSRLM] = {
                50, 5, 31, 7, 8, 4, 2.1, 0.5, 0.4, 2.9, 1.4
            };
            // Background uncertainties
            static const double BKGERR_LM[NUMSRLM] = {
                314, 22, 159, 37, 66, 33, 23, 3, 1, 1, 3
            };
            for (size_t ibin = 0; ibin < NUMSRLM; ++ibin)
            {
                stringstream ss; ss << "SR-LM-" << ibin+1;
                add_result(SignalRegionData(ss.str(), OBSNUM_LM[ibin], {_SRLM[ibin],  0.}, {BKGNUM_LM[ibin], BKGERR_LM[ibin]}));
            }

            // Observed event counts
            static const double OBSNUM_ML[NUMSRML] = {
                220, 2.5, 32, 0.9, 22, 9.5, 210, 36, 21.6, 11.6, 84, 15.5, 15.7, 5.3, 10.2, 27, 0.8, 17.8, 1.0, 17.8, 1.26, 830, 108, 117, 11.1, 111, 21, 42, 3.4, 320, 47, 51, 15.1, 131, 20, 27, 7.8, 12.9, 82, 18, 39, 4.9, 46, 5.7
            };
            // Background estimates
            static const double BKGNUM_ML[NUMSRML] = {
                40, 2.5, 6, 0.5, 4, 1.8, 40, 6, 3.2, 1.9, 11, 2.1, 2.2, 0.8, 2.1, 4, 0.5, 2.4, 0.4, 3.0, 0.33, 180, 22, 26, 2.3, 24, 5, 10, 0.9, 50, 8, 9, 2.6, 24, 4, 5, 1.5, 2.6, 14, 4, 8, 0.9, 10, 1.2
            };
            // Background uncertainties
            static const double BKGERR_ML[NUMSRML] = {
                263, 1, 34, 1, 28, 8, 265, 47, 20, 16, 105, 17, 21, 8, 12, 40, 2, 24, 0, 30, 2, 955, 136, 139, 8, 128, 20, 45, 3, 408, 50, 62, 24, 157, 24, 36, 11, 18, 117, 26, 29, 7, 44, 11
            };
            for (size_t ibin = 0; ibin < NUMSRML; ++ibin)
            {
                stringstream ss; ss << "SR-ML-" << ibin+1;
                add_result(SignalRegionData(ss.str(), OBSNUM_ML[ibin], {_SRML[ibin],  0.}, {BKGNUM_ML[ibin], BKGERR_ML[ibin]}));
            }
            
//            static const vector< vector<double> > BKGCOV = {
//                {},
//                {}
//            };

//            set_covariance(BKGCOV);

            return;
        }

    protected:
      void analysis_specific_reset() {
        for(size_t i=0;i<NUMSRHH;i++) { _SRHH[i]=0; }
        for(size_t i=0;i<NUMSRHL;i++) { _SRHL[i]=0; }
        for(size_t i=0;i<NUMSRLL;i++) { _SRLL[i]=0; }
        for(size_t i=0;i<NUMSRLM;i++) { _SRLM[i]=0; }
        for(size_t i=0;i<NUMSRML;i++) { _SRML[i]=0; }

      }

    };


    DEFINE_ANALYSIS_FACTORY(CMS_13TeV_2SSLEP_Stop_137invfb)

  }
}
