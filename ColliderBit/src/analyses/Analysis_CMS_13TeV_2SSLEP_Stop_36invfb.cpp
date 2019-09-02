///
///  \author Yang Zhang
///  \date 2019 July
///  *********************************************

// Based on http://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-16-035/index.html
//          https://arxiv.org/abs/1704.07323
            
// Search for physics beyond the standard model in events with two leptons of same sign, missing transverse momentum, and jets in proton-proton collisions at sqrt(s) = 13 TeV

// Note:
// 1. Not fully validated.
// 2. Not sure how to deal with the case that there are more than one lepton pairs. Here just use the first one.
// 3. Apply ATLAS Jet-Lepton Overlap Removal

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

    class Analysis_CMS_13TeV_2SSLEP_Stop_36invfb : public Analysis {
    protected:
        // Counters for the number of accepted events for each signal region
        static const size_t NUMSRHH = 51;
        double _SRHH[NUMSRHH];
        static const size_t NUMSRHL = 41;
        double _SRHL[NUMSRHL];
        static const size_t NUMSRLL = 8;
        double _SRLL[NUMSRLL];
        
        static const size_t NUMSRinc = 15;
        double _SRinc[NUMSRinc];
        static const size_t NUMSRexc = 15;
        double _SRexc[NUMSRexc];
        
        Cutflow _cutflow;
        
      // The following section copied from Analysis_ATLAS_1LEPStop_20invfb.cpp
      void JetLeptonOverlapRemoval(vector<HEPUtils::Jet*> &jetvec, vector<HEPUtils::Particle*> &lepvec, double DeltaRMax) {
        //Routine to do jet-lepton check
        //Discards jets if they are within DeltaRMax of a lepton

        vector<HEPUtils::Jet*> Survivors;

        for(unsigned int itjet = 0; itjet < jetvec.size(); itjet++) {
        bool overlap = false;
          HEPUtils::P4 jetmom=jetvec.at(itjet)->mom();
          for(unsigned int itlep = 0; itlep < lepvec.size(); itlep++) {
            HEPUtils::P4 lepmom=lepvec.at(itlep)->mom();
            double dR;

            dR=jetmom.deltaR_eta(lepmom);

            if(fabs(dR) <= DeltaRMax) overlap=true;
          }
          if(overlap) continue;
          Survivors.push_back(jetvec.at(itjet));
        }
        jetvec=Survivors;

        return;
      }

      void LeptonJetOverlapRemoval(vector<HEPUtils::Particle*> &lepvec, vector<HEPUtils::Jet*> &jetvec) {
        //Routine to do lepton-jet check
        //Discards leptons if they are within dR of a jet as defined in analysis paper

        vector<HEPUtils::Particle*> Survivors;

        for(unsigned int itlep = 0; itlep < lepvec.size(); itlep++) {
          bool overlap = false;
          HEPUtils::P4 lepmom=lepvec.at(itlep)->mom();
          for(unsigned int itjet= 0; itjet < jetvec.size(); itjet++) {
            HEPUtils::P4 jetmom=jetvec.at(itjet)->mom();
            double dR;
            double DeltaRMax = std::min(0.4, 0.04 + 10 / lepmom.pT());
            dR=jetmom.deltaR_eta(lepmom);

            if(fabs(dR) <= DeltaRMax) overlap=true;
          }
          if(overlap) continue;
          Survivors.push_back(lepvec.at(itlep));
        }
        lepvec=Survivors;

        return;
      }


    public:

        // Required detector sim
        static constexpr const char* detector = "CMS";

        Analysis_CMS_13TeV_2SSLEP_Stop_36invfb():
        _cutflow("CMS_13TeV_2SSLEP_Stop_36invfb", {"Trigger_and_2leptons", "At_least_one_SS_lepton_pair", "Baseline"})
        {

            set_analysis_name("CMS_13TeV_2SSLEP_Stop_36invfb");
            set_luminosity(36);
            
            for (size_t i = 0; i < NUMSRHH; ++i) _SRHH[i] = 0;
            for (size_t i = 0; i < NUMSRHL; ++i) _SRHL[i] = 0;
            for (size_t i = 0; i < NUMSRLL; ++i) _SRLL[i] = 0;
            
            for (size_t i = 0; i < NUMSRinc; ++i) _SRinc[i] = 0;
            for (size_t i = 0; i < NUMSRexc; ++i) _SRexc[i] = 0;
            
            
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
            const vector<double> bMu={0., 10, 15., 20., 25, 30, 40, 50, DBL_MAX};  // Bin edges in pT. Assume flat efficiency above 200, where the CMS map stops.
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

            // Overlap removal
            JetLeptonOverlapRemoval(candJets,electrons,0.2);
            LeptonJetOverlapRemoval(electrons,candJets);
            JetLeptonOverlapRemoval(candJets,muons,0.4);
            LeptonJetOverlapRemoval(muons,candJets);

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
                else if( jet->pT() > 40. ) {
                    nonbJets.push_back(jet);
                }
            }
            
            size_t Nb=bJets.size();
            size_t Nj=nonbJets.size();
            
            // Leptons = electrons + muons
            vector<HEPUtils::Particle*> leptons;
            leptons=electrons;
            leptons.insert(leptons.end(),muons.begin(),muons.end());
            sort(leptons.begin(),leptons.end(),comparePt);

            // At least two light leptons
            if (leptons.size()<2) return;
            
            // Triggers
            bool pure_dilepton_trigger=false;
            // Leading electron (muon) PT > 23 (17) GeV
            // Subleading electron (muon) PT > 12 (8) GeV
            if (leptons[0]->pT() > ( abs(leptons[0]->pid()) == 11 ? 23 : 17 ) \
            and leptons[1]->pT() > ( abs(leptons[1]->pid()) == 11 ? 12 : 8 ) ){
                pure_dilepton_trigger = true;
            }
            if ( not pure_dilepton_trigger and HT<300 ) return;
            _cutflow.fill(1); // Trigger and >=2 leptons
            
            // Find pair same sign (SS) leptons
            vector<size_t> SS_1,SS_2;
            for (size_t i=0; i<leptons.size(); ++i) {
                for (size_t j=i+1; j<leptons.size(); ++j) {
                    if (leptons[i]->pid()*leptons[j]->pid()>0 and (leptons[i]->mom()+leptons[j]->mom()).m()>8){
                        SS_1.push_back(i);
                        SS_2.push_back(j);
                    }
                }
            }

            // At least one SS lepton pair ( with an invari-ant mass above 8 GeV )
            if (SS_1.size()==0) return;
            _cutflow.fill(2); // At least one SS lepton pair
            
            // An additional loose lepton forms an opposite-sign same-flavor pair
            // withone of the two SS leptons, with an invariant mass less than 12 GeV 
            // or between 76 and 106 GeV
            if (leptons.size()>2){
                for (size_t i=0; i<SS_1.size(); ++i) {
                    for (size_t j=0; j<leptons.size(); ++j) {
                        if ( j != SS_1[i] and j != SS_2[i]) {
                            if (leptons[j]->pid()+leptons[SS_1[i]]->pid()==0){
                                double mll_additional = (leptons[j]->mom()+leptons[SS_1[i]]->mom()).m();
                                if ( mll_additional < 12 or (mll_additional>76 and mll_additional<106)) return;
                            }
                            if (leptons[j]->pid()+leptons[SS_2[i]]->pid()==0){
                                double mll_additional = (leptons[j]->mom()+leptons[SS_2[i]]->mom()).m();
                                if ( mll_additional < 12 or (mll_additional>76 and mll_additional<106)) return;
                            }
                        }
                    }
                }
            }
            
            
            // At least two jets and MET>50
            if ( nonbJets.size()<2 or  met<50) return;
            _cutflow.fill(3); // Baseline (two jets and MET>50 GeV)
            
            // M_T^{miss}
            // The smallest of the transverse masses constructed between p^miss_T and each of the leptons.
            double MTmiss = 9999;
            for (HEPUtils::Particle* lep : leptons) {
                double MTmiss_temp = sqrt(2.*lep->pT()*met*(1. - cos(lep->mom().deltaPhi(ptot))));
                if (MTmiss_temp<MTmiss) {
                    MTmiss = MTmiss_temp;
                }
            }

            bool pp = leptons[SS_1[0]]->pid()>0; // TODO: Not sure which lepton pair.
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
            bool SSHH_combine = (MTmiss_l_120 and met_50_200   and Nj_5 and HT_300) or \
                                (MTmiss_l_120 and met_200_300           and HT_300) or \
                                (MTmiss_g_120 and met<300               and HT_300);

            // SR HH
            if ( leptons[SS_1[0]]->pT() > 25. and leptons[SS_2[0]]->pT() > 25.) {
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
                    if (MTmiss_l_120 and met_50_200             and HT_300_1125 and pp)  _SRHH[36]++;
                    if (MTmiss_l_120 and met_50_200             and HT_300_1125 and !pp) _SRHH[37]++;
                    if (MTmiss_l_120 and met_200_300            and HT_300_1125)         _SRHH[38]++;
                    if (MTmiss_g_120 and met<300                and HT_300)              _SRHH[39]++;
                    if (MTmiss_g_120 and met<300                and HT_300_1125)         _SRHH[40]++;   
                }
                
                if (met_300_500 and HT>300       and pp)  _SRHH[41]++;
                if (met_300_500 and HT>300       and !pp) _SRHH[42]++;
                if (met_500     and HT>300       and pp)  _SRHH[43]++;
                if (met_500     and HT>300       and !pp) _SRHH[44]++;
                
                if (met<300     and HT_1125_1300 and pp)  _SRHH[45]++;
                if (met<300     and HT_1125_1300 and !pp) _SRHH[46]++;
                if (met<300     and HT_1300_1600 and pp)  _SRHH[47]++;
                if (met<300     and HT_1300_1600 and !pp) _SRHH[48]++;
                if (met<300     and HT_1600 and pp)  _SRHH[48]++;
                if (met<300     and HT_1600 and !pp) _SRHH[50]++;
                
            }
            
            bool SSHL_combine = MTmiss_l_120&&( (met_50_200&&Nj_5) or met_200_300 )&&HT_300 ;
            
            // SR HL
            if ( leptons[SS_1[0]]->pT() > 25. and leptons[SS_2[0]]->pT() < 25.) {
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
                    if ( met_200_300 and Nj_2_4 and HT_300_1125 and pp)  _SRHL[13]++;
                    if ( met_200_300 and Nj_2_4 and HT_300_1125 and !pp) _SRHL[14]++;
                    
                    if ( met_200_300 and Nj_5 and HT_300_1125 and pp)    _SRHL[15]++;
                    if ( met_200_300 and Nj_5 and HT_300_1125 and !pp)   _SRHL[16]++;
                } else if(Nb==2 and MTmiss_l_120) {
                    if ( met_50_200  and Nj_2_4 and HT_300)              _SRHL[17]++;
                    if ( met_50_200  and Nj_2_4 and HT_300_1125)         _SRHL[18]++;
                    if ( SSHL_combine and pp)                            _SRHL[19]++;
                    if ( SSHL_combine and !pp)                           _SRHL[20]++;
                    if ( met_50_200  and Nj_5 and HT_300_1125 and pp)    _SRHL[21]++;
                    if ( met_50_200  and Nj_5 and HT_300_1125 and !pp)   _SRHL[22]++;
                    if ( met_200_300 and Nj_2_4 and HT_300_1125 and pp)  _SRHL[23]++;
                    if ( met_200_300 and Nj_2_4 and HT_300_1125 and !pp) _SRHL[24]++;
                    if ( met_200_300 and Nj_5 and HT_300_1125)           _SRHL[25]++;
                }else if(Nb==3 and MTmiss_l_120) {
                    if ( met_50_200 and HT_300 and pp)                   _SRHL[26]++;
                    if ( met_50_200 and HT_300 and !pp)                  _SRHL[27]++;
                    if ( met_50_200 and HT_300_1125 and pp)              _SRHL[28]++;
                    if ( met_50_200 and HT_300_1125 and !pp)             _SRHL[29]++;
                    if ( met_200_300 and HT_300_1125)                    _SRHL[30]++;
                }
                if (MTmiss_g_120 and met<300 and HT_300)       _SRHL[31]++;
                if (MTmiss_g_120 and met<300 and HT_300_1125)  _SRHL[32]++;
                
                if (met_300_500  and HT>300 and pp)  _SRHL[33]++;
                if (met_300_500  and HT>300 and !pp) _SRHL[34]++;
                if (met_500      and HT>300 and pp)  _SRHL[35]++;
                if (met_500      and HT>300 and !pp) _SRHL[36]++;
                
                if (met<300      and HT_1125_1300 and pp)  _SRHL[37]++;
                if (met<300      and HT_1125_1300 and !pp) _SRHL[38]++;
                if (met<300      and HT_1300 and pp)       _SRHL[39]++;
                if (met<300      and HT_1300 and !pp)      _SRHL[40]++;
            }

            // SR LL
            if (leptons[SS_1[0]]->pT() < 25. and leptons[SS_2[0]]->pT() < 25.) {
                if (HT>300) {
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
            
            // Inclusive SR
            if (  leptons[SS_1[0]]->pT() > 25. and leptons[SS_2[0]]->pT() > 25. ) {
                // Nj>=2 and met>50 have been applied
                if ( Nb==0 and HT>=1200)               _SRinc[0]++;
                if ( Nb>=2 and HT>=1100)               _SRinc[1]++;
                if ( Nb==0 and met>450)                _SRinc[2]++;
                if ( Nb>=2 and met>300)                _SRinc[3]++;
                if ( Nb==0 and met>250 and MTmiss>120) _SRinc[4]++;
                if ( Nb>=2 and met>150 and MTmiss>120) _SRinc[5]++;
                if ( Nb==0 and HT>900 and met>200)     _SRinc[6]++;
                if ( Nb>=2 and HT>900 and met>200)     _SRinc[7]++;
                if ( Nj>=7)                            _SRinc[8]++;
                if ( Nj>=4 and MTmiss>120)             _SRinc[9]++;
                if ( Nb>=3)                            _SRinc[10]++;
                if ( HT>700)                           _SRinc[11]++;
            }
            
            if (  leptons[SS_1[0]]->pT() < 25. and leptons[SS_2[0]]->pT() < 25. ) {
                // Nj>=2 and met>50 have been applied
                if (met>200) _SRinc[12]++;
                if (Nj>=5)   _SRinc[13]++;
                if (Nb>=3)   _SRinc[14]++;
            }
            
            // Exclusive SR
            if (  leptons[SS_1[0]]->pT() > 25. and leptons[SS_2[0]]->pT() > 25. ) {
                // Nj>=2 and met>50 have been applied
                if (Nb==0 and met<300 and HT<1125 and (HT<300 or MTmiss<120)) _SRexc[0]++;
                if (Nb==0 and met<300 and HT<1125 and HT>300 and MTmiss>120)  _SRexc[1]++;
                if (Nb==1 and met<300 and HT<1125 and (HT<300 or MTmiss<120)) _SRexc[2]++;
                if (Nb==1 and met<300 and HT<1125 and HT>300 and MTmiss>120)  _SRexc[3]++;
                if (Nb==2 and met<300 and HT<1125 and (HT<300 or MTmiss<120)) _SRexc[4]++;
                if (Nb==2 and met<300 and HT<1125 and HT>300 and MTmiss>120)  _SRexc[5]++;
                if (Nb>=3 and met<300 and HT<1125 and (HT<300 or MTmiss<120)) _SRexc[6]++;
                if (Nb>=3 and met<300 and HT<1125 and HT>300 and MTmiss>120)  _SRexc[7]++;
                if (          met>300 and             HT>300)                 _SRexc[8]++;
                if (          met<300 and HT>1125)                            _SRexc[9]++;
            }
            
            if (  leptons[SS_1[0]]->pT() > 25. and leptons[SS_2[0]]->pT() < 25. ) {
                // Nj>=2 and met>50 have been applied
                if (met<300 and HT<1125 and MTmiss<120) _SRexc[10]++;
                if (met<300 and HT<1125 and MTmiss>120) _SRexc[11]++;
                if (met>300 and HT>300)                 _SRexc[12]++;
                if (met<300 and HT>1125)                _SRexc[13]++;
            }
            if (  leptons[SS_1[0]]->pT() < 25. and leptons[SS_2[0]]->pT() < 25. ) {
                // Nj>=2 and met>50 have been applied
                if (HT>300) _SRexc[14]++;
            }
            
            return;
        }

        /// Combine the variables of another copy of this analysis (typically on another thread) into this one.
        void combine(const Analysis* other)
        {
            const Analysis_CMS_13TeV_2SSLEP_Stop_36invfb* specificOther
                = dynamic_cast<const Analysis_CMS_13TeV_2SSLEP_Stop_36invfb*>(other);
            
            for (size_t i = 0; i < NUMSRHH; ++i) _SRHH[i] += specificOther->_SRHH[i];
            for (size_t i = 0; i < NUMSRHL; ++i) _SRHL[i] += specificOther->_SRHL[i];
            for (size_t i = 0; i < NUMSRLL; ++i) _SRLL[i] += specificOther->_SRLL[i];
            
            for (size_t i = 0; i < NUMSRinc; ++i) _SRinc[i] += specificOther->_SRinc[i];
            for (size_t i = 0; i < NUMSRexc; ++i) _SRexc[i] += specificOther->_SRexc[i];
            
        }


        void collect_results() {

            #ifdef CHECK_CUTFLOW
            cout << _cutflow << endl;
            #endif

            // Observed event counts
            static const double OBSNUM_HH[NUMSRHH] = {
                468, 162, 24.4, 17.6, 17.8, 7.8, 1.96, 4.58, 3.63, 2.82, 313, 104, 9.5, 8.7, 14.4, 12.7, 7.3, 3.92, 3.26, 2.6, 3.02, 2.80, 70, 35.7, 3.99, 2.68, 9.7, 7.9, 2.78, 1.86, 2.20, 1.85, 1.20, 1.81, 1.98, 1.43, 4.2, 3.04, 0.63, 0.29, 0.80, 13.4, 8.0, 3.33, 0.94, 2.92, 1.78, 1.95, 1.23, 1.46, 0.74
            };
            // Background estimates
            static const double BKGNUM_HH[NUMSRHH] = {
                98, 25, 5.4, 3.0, 3.9, 1.5, 0.47, 0.81, 0.75, 0.56, 87, 20, 1.9, 2.0, 2.9, 2.6, 1.2, 0.79, 0.74, 2.7, 0.75, 0.57, 12, 5.9, 0.73, 0.80, 1.8, 2.5, 0.58, 0.38, 0.54, 0.39, 0.32, 0.42, 0.61, 0.37, 1.3, 0.68, 0.17, 0.34, 0.22, 1.9, 3.0, 0.74, 0.26, 0.50, 0.42, 0.39, 0.30, 0.31, 0.18
            };
            // Background uncertainties
            static const double BKGERR_HH[NUMSRHH] = {
                435, 166, 30, 24, 22, 6, 2, 5, 3, 3, 304, 111, 13, 11, 17, 10, 11, 2, 3, 4, 3, 1, 90, 40, 2, 0, 9, 8, 1, 1, 1, 5, 0, 3, 1, 2, 2, 4, 1, 0, 3, 19, 8, 3, 1, 3, 3, 5, 3, 0, 0
            };
            for (size_t ibin = 0; ibin < NUMSRHH; ++ibin)
            {
                stringstream ss; ss << "SR-HH-" << ibin+1;
                add_result(SignalRegionData(ss.str(), OBSNUM_HH[ibin], {_SRHH[ibin],  0.}, {BKGNUM_HH[ibin], BKGERR_HH[ibin]}));
                #ifdef CHECK_CUTFLOW
                cout << ss.str() << "\t" << _SRHH[ibin] << endl;
                #endif
            }

            // Observed event counts
            static const double OBSNUM_HL[NUMSRHL] = {
                419, 100, 9.2, 15.0, 7.3, 4.1, 1.01, 300, 73, 2.30, 2.24, 12.8, 8.9, 4.5, 4.7, 2.3, 0.73, 54, 23.7, 0.59, 0.34, 5.2, 4.9, 0.97, 1.79, 1.01, 1.03, 1.33, 2.89, 2.24, 0.27, 0.79, 0.53, 6.3, 2.92, 0.51, 0.15, 1.07, 0.81, 1.54, 1.23
            };
            // Background estimates
            static const double BKGNUM_HL[NUMSRHL] = {
                100, 20, 2.4, 4.5, 1.5, 1.2, 0.28, 82, 17, 0.61, 0.87, 3.3, 2.3, 1.3, 1.6, 1.1, 0.29, 12, 4.9, 0.17, 0.20, 1.2, 1.4, 0.27, 0.74, 0.27, 0.44, 0.61, 0.99, 0.79, 0.30, 0.33, 0.13, 1.3, 0.87, 0.15, 0.07, 0.33, 0.47, 0.50, 0.53
            };
            // Background uncertainties
            static const double BKGERR_HL[NUMSRHL] = {
                442, 101, 6, 13, 14, 5, 0, 346, 95, 1, 1, 12, 8, 5, 4, 1, 1, 62, 24, 2, 1, 9, 6, 0, 0, 1, 1, 0, 3, 2, 1, 1, 0, 6, 3, 3, 0, 3, 0, 4, 1
            };
            for (size_t ibin = 0; ibin < NUMSRHL; ++ibin)
            {
                stringstream ss; ss << "SR-HL-" << ibin+1;
                add_result(SignalRegionData(ss.str(), OBSNUM_HL[ibin], {_SRHL[ibin],  0.}, {BKGNUM_HL[ibin], BKGERR_HL[ibin]}));
                #ifdef CHECK_CUTFLOW
                cout << ss.str() << "\t" << _SRHL[ibin] << endl;
                #endif
            
            }

            // Observed event counts
            static const double OBSNUM_LL[NUMSRLL] = {
                12.0, 1.88, 15.5, 1.42, 4.2, 0.84, 0.95, 0.09
            };
            // Background estimates
            static const double BKGNUM_LL[NUMSRLL] = {
                3.9, 0.62, 4.7, 0.69, 1.4, 0.48, 0.52, 0.07
            };
            // Background uncertainties
            static const double BKGERR_LL[NUMSRLL] = {
                12, 3, 17, 4, 5, 2, 0, 0
            };
            for (size_t ibin = 0; ibin < NUMSRLL; ++ibin)
            {
                stringstream ss; ss << "SR-LL-" << ibin+1;
                add_result(SignalRegionData(ss.str(), OBSNUM_LL[ibin], {_SRLL[ibin],  0.}, {BKGNUM_LL[ibin], BKGERR_LL[ibin]}));
                #ifdef CHECK_CUTFLOW
                cout << ss.str() << "\t" << _SRLL[ibin] << endl;
                #endif
            }

            return;
        }

    protected:
      void analysis_specific_reset() {
        for(size_t i=0;i<NUMSRHH;i++) { _SRHH[i]=0; }
        for(size_t i=0;i<NUMSRHL;i++) { _SRHL[i]=0; }
        for(size_t i=0;i<NUMSRLL;i++) { _SRLL[i]=0; }
        
        for(size_t i=0;i<NUMSRinc;i++) { _SRinc[i]=0; }
        for(size_t i=0;i<NUMSRexc;i++) { _SRexc[i]=0; }
        
      }

    };


    DEFINE_ANALYSIS_FACTORY(CMS_13TeV_2SSLEP_Stop_36invfb)

    //
    // Derived analysis class for the 2Lep0Jets SRs
    //
    class Analysis_CMS_13TeV_2SSLEP_Stop_inclusive_36invfb : public Analysis_CMS_13TeV_2SSLEP_Stop_36invfb {

    public:
      Analysis_CMS_13TeV_2SSLEP_Stop_inclusive_36invfb() {
        set_analysis_name("CMS_13TeV_2SSLEP_Stop_inclusive_36invfb");
      }

        virtual void collect_results() {
                  
            // Observed event counts
            static const double OBSNUM_inc[NUMSRinc] = {
                4.00, 3.63, 3.72, 3.32, 1.68, 3.82, 5.6, 5.8, 10.1, 15.2, 13.3, 3.6, 4.9, 7.3, 1.06
            };
            // Background estimates
            static const double BKGNUM_inc[NUMSRinc] = {
                0.79, 0.71, 0.83, 0.81, 0.44, 0.76, 1.1, 1.3, 2.7, 3.5, 3.4, 2.5, 2.9, 5.5, 0.99
            };
            // Background uncertainties
            static const double BKGERR_inc[NUMSRinc] = {
                10, 4, 4, 6, 2, 7, 10, 9, 9, 22, 17, 3, 10, 6, 0
            };
            for (size_t ibin = 0; ibin < NUMSRinc; ++ibin)
            {
                stringstream ss; ss << "InSR-" << ibin+1;
                add_result(SignalRegionData(ss.str(), OBSNUM_inc[ibin], {_SRinc[ibin],  0.}, {BKGNUM_inc[ibin], BKGERR_inc[ibin]}));
                #ifdef CHECK_CUTFLOW
                cout << ss.str() << "\t" << _SRinc[ibin] << endl;
                #endif
            
            }
      
        }

    };

    // Factory fn
    DEFINE_ANALYSIS_FACTORY(CMS_13TeV_2SSLEP_Stop_inclusive_36invfb)

    //
    // Derived analysis class for the 2Lep0Jets SRs
    //
    class Analysis_CMS_13TeV_2SSLEP_Stop_exclusive_36invfb : public Analysis_CMS_13TeV_2SSLEP_Stop_36invfb {

    public:
      Analysis_CMS_13TeV_2SSLEP_Stop_exclusive_36invfb() {
        set_analysis_name("CMS_13TeV_2SSLEP_Stop_exclusive_36invfb");
      }

        virtual void collect_results() {
                  
            // Observed event counts
            static const double OBSNUM_exc[NUMSRexc] = {
                700, 11.0, 477, 8.4, 137, 4.9, 11.6, 0.8, 25.7, 10.1, 1070, 1.33, 9.9, 4.7, 37
            };
            // Background estimates
            static const double BKGNUM_exc[NUMSRexc] = {
                130, 2.2, 120, 3.5, 25, 1.2, 3.1, 0.24, 5.4, 2.2, 250, 0.46, 2.5, 1.8, 12
            };
            // Background uncertainties
            static const double BKGERR_exc[NUMSRexc] = {
                685, 11, 482, 8, 152, 8, 10, 3, 31, 14, 1167, 1, 12, 8, 43
            };
            for (size_t ibin = 0; ibin < NUMSRexc; ++ibin)
            {
                stringstream ss; ss << "ExSR-" << ibin+1;
                add_result(SignalRegionData(ss.str(), OBSNUM_exc[ibin], {_SRexc[ibin],  0.}, {BKGNUM_exc[ibin], BKGERR_exc[ibin]}));
                #ifdef CHECK_CUTFLOW
                cout << ss.str() << "\t" << _SRexc[ibin] << endl;
                #endif
            
            }

            static const vector< vector<double> > BKGCOV = {
                {17559.1, 111.4, 13059.1, 136.0, 1982.8, 34.0, 187.1, 14.8, 260.9, 80.2, 26183.7, 21.0, 134.9, 78.2, 987.5},
                {111.4, 3.8, 85.1, 1.9, 20.0, 0.6, 2.5, 0.2, 3.9, 1.7, 186.6, 0.2, 1.3, 0.8, 6.6},
                {13059.1, 85.1, 12489.9, 102.2, 1847.0, 32.4, 158.8, 11.1, 179.7, 58.6, 23663.6, 19.6, 123.0, 63.8, 881.1},
                {136.0, 1.9, 102.2, 10.4, 28.8, 1.1, 3.8, 0.3, 6.3, 2.3, 210.4, 0.3, 1.9, 1.4, 8.5},
                {1982.8, 20.0, 1847.0, 28.8, 525.1, 11.2, 36.6, 2.5, 43.8, 18.3, 3573.7, 3.7, 24.4, 12.8, 139.8},
                {34.0, 0.6, 32.4, 1.1, 11.2, 1.1, 1.3, 0.1, 1.8, 0.7, 60.2, 0.1, 0.7, 0.4, 2.4},
                {187.1, 2.5, 158.8, 3.8, 36.6, 1.3, 9.1, 0.3, 6.3, 2.4, 320.9, 0.3, 2.3, 1.4, 12.5},
                {14.8, 0.2, 11.1, 0.3, 2.5, 0.1, 0.3, 0.1, 0.5, 0.2, 21.4, 0.0, 0.2, 0.1, 0.8},
                {260.9, 3.9, 179.7, 6.3, 43.8, 1.8, 6.3, 0.5, 25.7, 3.6, 345.5, 0.4, 2.9, 2.3, 14.5},
                {80.2, 1.7, 58.6, 2.3, 18.3, 0.7, 2.4, 0.2, 3.6, 4.0, 126.8, 0.2, 1.2, 0.9, 4.2},
                {26183.7, 186.6, 23663.6, 210.4, 3573.7, 60.2, 320.9, 21.4, 345.5, 126.8, 60200.4, 45.8, 307.4, 170.1, 1998.6},
                {21.0, 0.2, 19.6, 0.3, 3.7, 0.1, 0.3, 0.0, 0.4, 0.2, 45.8, 0.2, 0.3, 0.1, 1.7},
                {134.9, 1.3, 123.0, 1.9, 24.4, 0.7, 2.3, 0.2, 2.9, 1.2, 307.4, 0.3, 6.2, 1.1, 11.2},
                {78.2, 0.8, 63.8, 1.4, 12.8, 0.4, 1.4, 0.1, 2.3, 0.9, 170.1, 0.1, 1.1, 3.3, 6.0},
                {987.5, 6.6, 881.1, 8.5, 139.8, 2.4, 12.5, 0.8, 14.5, 4.2, 1998.6, 1.7, 11.2, 6.0, 135.5}
            };

            set_covariance(BKGCOV);
      
        }

    };

    // Factory fn
    DEFINE_ANALYSIS_FACTORY(CMS_13TeV_2SSLEP_Stop_exclusive_36invfb)

  }
}
