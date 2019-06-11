///
///  \author Tomas Gonzalo
///  \date 2019 June
///  *********************************************

// Based on https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/CONFNOTES/ATLAS-CONF-2019-016/

// - 139 fb^-1 data

#include <vector>
#include <cmath>
#include <memory>
#include <iomanip>
#include <algorithm>
#include <fstream>

#include "gambit/ColliderBit/analyses/Analysis.hpp"
#include "gambit/ColliderBit/ATLASEfficiencies.hpp"
#include "gambit/ColliderBit/mt2_bisect.h"

// #define CHECK_CUTFLOW

using namespace std;

namespace Gambit
{
  namespace ColliderBit
  {

    class Analysis_ATLAS_13TeV_2OSLEP_Z_139invfb : public Analysis
    {

    protected:

      // Counters for the number of accepted events for each signal region
      std::map<string,double> _numSR = {
        {"SR1A", 0},
        {"SR1B", 0},
        {"SR2A", 0},
        {"SR2B", 0},
      };

    private:

      struct ptComparison
      {
        bool operator() (HEPUtils::Particle* i,HEPUtils::Particle* j) {return (i->pT()>j->pT());}
      } comparePt;

      struct ptJetComparison
      {
        bool operator() (HEPUtils::Jet* i,HEPUtils::Jet* j) {return (i->pT()>j->pT());}
      } compareJetPt;

      // Jet lepton overlap removal
      // Discards jets if they are within DeltaRMax of a lepton
      void JetLeptonOverlapRemoval(vector<HEPUtils::Jet*>& jets, vector<HEPUtils::Particle*>& leptons, double DeltaRMax, double pTMax = 1000., double btaggeff = 1)
      {
        vector<HEPUtils::Jet*> jetSurvivors;
        vector<HEPUtils::Particle*> leptonSurvivors;
        for(HEPUtils::Jet* jet : jets)
        {
          bool overlap = false;
          for(HEPUtils::Particle* lepton : leptons)
          {
            double dR = jet->mom().deltaR_eta(lepton->mom());
            if(jet->pT() < pTMax && fabs(dR) <= DeltaRMax) overlap = true;
            // TODO: Add conditional for b-tagging efficiency
          }
          if(!overlap) jetSsurvivors.push_back(jet);
        }
        jets = jetSurvivors;
        return;
      }

      // Lepton jet overlap removal
      // Discards leptons if they are within DeltaRMax of a jet
      void LeptonJetOverlapRemoval(vector<HEPUtils::Particle*>& leptons, vector<HEPUtils::Jet*>& jets, double DeltaRMax)
      {
        vector<HEPUtils::Particle*> survivors;
        for(HEPUtils::Particle* lepton : leptons)
        {
          bool overlap = false;
          for(HEPUtils::Jet* jet : jets)
          {
            double dR = jet->mom().deltaR_eta(lepton->mom());
            if(fabs(dR) <= DeltaRMax) overlap = true;
          }
          if(!overlap) survivors.push_back(lepton);
        }
        leptons = survivors;
        return;
      }

      // Particle overlap removal
      // Discards particle (from "particles1") if it is within DeltaRMax of another particle
      void ParticleOverlapRemoval(vector<HEPUtils::Particle*>& particles1, vector<HEPUtils::Particle*>& particles2, double DeltaRMax)
      {
        vector<HEPUtils::Particle*> survivors;
        for(HEPUtils::Particle* p1 : particles1)
        {
          bool overlap = false;
          for(HEPUtils::Particle* p2 : particles2)
          {
            double dR = p1->mom().deltaR_eta(p2->mom());
            if(fabs(dR) <= DeltaRMax) overlap = true;
          }
          if(!overlap) survivors.push_back(p1);
        }
        particles1 = survivors;
        return;
      }

      // Removes a lepton from the leptons1 vector if it forms an OS pair with a
      // lepton in leptons2 and the pair has a mass in the range (m_low, m_high).
      void removeOSPairsInMassRange(vector<HEPUtils::Particle*>& leptons1, vector<HEPUtils::Particle*>& leptons2, double m_low, double m_high)
      {
        vector<HEPUtils::Particle*> l1_survivors;
        for(HEPUtils::Particle* l1 : leptons1)
        {
          bool survived = true;
          for(HEPUtils::Particle* l2 : leptons2)
          {
            if(l2 == l1) continue;
            if (l1->pid()*l2->pid() < 0.)
            {
              double m = (l1->mom() + l2->mom()).m();
              if ((m >= m_low) && (m <= m_high))
              {
                survived = false;
                break;
              }
            }
          }
          if(survived) l1_survivors.push_back(l1);
        }
        leptons1 = l1_survivors;
        return;
      }

    public:

      // Required detector sim
      static constexpr const char* detector = "ATLAS";

      Analysis_ATLAS_13TeV_2OSLEP_Z_139invfb()
      {

        set_analysis_name("ATLAS_13TeV_2OSLEP_Z_139invfb");
        set_luminosity(139);

        #ifdef CHECK_CUTFLOW
          NCUTS = 11;
          for (size_t i=0;i<NCUTS;i++)
          {
            cutFlowVector.push_back(0);
            cutFlowVectorATLAS_400_0.push_back(0);
            cutFlowVector_str.push_back("");
          }
        #endif

      }

      void run(const HEPUtils::Event* event)
      {

        // Baseline objects
        vector<HEPUtils::Particle*> baselineElectrons;
        vector<HEPUtils::Particle*> baselineMuons;
        vector<HEPUtils::Particle*> baselineTaus;
        vector<HEPUtils::Jet*> baselineJets;
        double met = event->met();

        // Electron candidates are reconstructed from isolated electromagnetic calorimeter energy deposits matched to ID tracks and are required to have |η| < 2.47, a transverse momentum pT > 4.5 GeV, and to pass the “LooseAndBLayer” requirement in arXiv: 1902.04655 [hep-ex].
       for (HEPUtils::Particle* electron : event->electrons())
        {
          if (electron->pT()>4.5 && electron->abseta()<2.47) baselineElectrons.push_back(electron);
        }

        // Apply electron efficiency
        ATLAS::applyElectronEff(baselineElectrons);

        // Apply loose electron selection
        // TODO: Check that the LooseAndBLayer cut is the same as this
        ATLAS::applyLooseIDElectronSelectionR2(baselineElectrons);

        // Muon candidates are reconstructed in the region |η| < 2.4 from muon spectrometer tracks matching ID tracks. Candidate muons must have pT > 4 GeV and pass the medium identification requirements defined in arXiv: 1603.05598 [hep-ex]. 
        for (HEPUtils::Particle* muon : event->muons())
        {
          if (muon->pT()>4. && muon->abseta()<2.4) baselineMuons.push_back(muon);
        }

        // Apply muon efficiency
        ATLAS::applyMuonEff(baselineMuons);

        // TODO Apply "medium" muon ID criteria

        // TODO: transverse and longitudinal impact parameter cuts

        // Only jet candidates with pT > 20 GeV and |η| < 2.8 are considered in the analysis
        for (HEPUtils::Jet* jet : event->jets())
        {
          if (jet->pT()>20. && jet->abseta()<2.8) baselineJets.push_back(jet);
        }
        // TODO: Additional requirements for jets

        // TODO: Something about b-jets

        // Overlap removal

        // 1) Remove jets within DeltaR = 0.2 of electron
        // TODO: remove electron if b-taggin effeciciency > 85%
        JetLeptonOverlapRemoval(baselineJets, baselineElectrons, 0.2, 200., 0.85);

        // 2) Remove electrons within DeltaR = 0.4 of a jet
        LeptonJetOverlapRemoval(baselineElectrons, baselineJets, 0.4);

        // 3) Remove jets within DeltaR = 0.2 of a muon
        JetLeptonOverlapRemoval(baselineJets, baselineMuons, 0.2);

        // 4) Remove muons within DeltaR = 0.4 of jet
        // TODO: Not quite 0.4 but min(0.4, 0.04 + pT(µ)/10 GeV)
        // TODO: Remove the jet instead if the jet has fewer than 3 associated tracks
        LeptonJetOverlapRemoval(baselineMuons, baselineJets, 0.4);

        // 5) Remove electron candidates sharing and ID track with a muon candidate
        // TODO: Missing


        // Signal objects
        vector<HEPUtils::Jet*> signalJets = baselineJets;
        vector<HEPUtils::Particle*> signalElectrons = baselineElectrons;
        vector<HEPUtils::Particle*> signalMuons;
        vector<HEPUtils::Particle*> signalLeptons;

         // TODO: Signal electrons must satisfy the “medium” identification requirement as defined in arXiv: 1902.04655 [hep-ex]

        // Signal muons must have pT > 5 GeV.
        for (HEPUtils::Particle* signalMuon : baselineMuons)
        {
          if (signalMuon->pT() > 5.) signalMuons.push_back(signalMuon);
        }
         
        // TODO: isolation criteria for signal leptons (see paper)

        // TODO: Check that nothing to be done for corrections between MC and signal

        // Fill signal leptons
        signalLeptons = signalElectrons;
        signalLeptons.insert(signalLeptons.end(), signalMuons.begin(), signalMuons.end());

        // Sort by pT
        sort(signalJets.begin(), signalJets.end(), compareJetPt);
        sort(signalLeptons.begin(), signalLeptons.end(), comparePt);

        // Trigger requirements are
        // - >=3 signal leptons
        // - >=1 SF-OS pair
        // - leading lepton pT > 40 GeV
        // - subleading lepton pT > 20 GeV
        // - Zlike, |mll - mZ| < 15 GeV
        bool trigger = false;

        // Count signal leptons and jets
        size_t nSignalLeptons = signalLeptons.size();
        size_t nSignalJets = signalJets.size();

        // Get SFOS pairs
        vector<vector<HEPUtils::Particle*>> SFOSpairs = getSFOSpairs(signalLeptons);

        // Get SFOS pairs masses and pTs
        vector<double> SFOSpair_masses;
        vector<double> SFOSpair_pTs;
        for (vector<HEPUtils::Particle*> pair : SFOSpairs)
        {
          SFOSpair_masses.push_back( (pair.at(0)->mom() + pair.at(1)->mom()).m() );
          SFOSpair_pTs.push_back( (pair.at(0)->mom() + pair.at(1)->mom()).pT() );
        }
        std::sort(SFOSpair_masses.begin(), SFOSpair_masses.end(), std::greater<double>());
        std::sort(SFOSpair_pTs.begin(), SFOSpair_pTs.end(), std::greater<double>());

        // Z resonance
        bool Zlike = false;
        double mZ = 91.2;
        for(double m : SFOSpair_masses)
        {
          if (abs(m - mZ) < 15)
            Zlike = true;
        }

        // Combine all trigger cuts
        trigger = nSignalLeptons >= 3 && SFOSpairs.size() >= 1 && signalLeptons.at(0)->pT() > 40. && signalLeptons.at(1)->pT() > 20. && Zlike;     


        // Effective mass (met + pT of all signal leptons + pT of all jets with pT>40 GeV)
//        double meff = met;
//        for (HEPUtils::Particle* l : signalLeptons)
//        {
//          meff += l->pT();
//        }
//        for (HEPUtils::Jet* jet : signalJets)
//        {
//          if(jet->pT()>40.) meff += jet->pT();
//        }


        // Signal Regions

        // Requirement                      SR1A     SR1B    SR2A    SR2B
        // ---------------------------------------------------------------
        // Third leading lepton pT           >20      >20     <20     <60   // done
        // njets (pT > 30 GeV)               >=4      >=5     >=3     >=3   // done
        // nb-tagged jets (pT > 30 GeV)      >=1      >=1      -      >=1   // TODO
        // Leading jet pT                     -        -     >150      -    // done
        // Leading b-tagged jet pT            -      >100      -       -    // TODO
        // MET                              >250     >150    >200    >350   // done
        // pTll                               -      >150     <50    >150   // done
        // mT23l                            >100       -       -       -    // TODO

        // SR1A
        if (trigger && 
            signalLeptons.at(2)->pT() > 20. && 
            nSignalJets >= 4 && signalJets.at(3)->pT() > 30. &&
            // TODO
            // -
            // TODO
            met > 250. &&
            // -
            // TODO
           ) 
          _numSR["SR1A"]++;

        // SR1B
        if (trigger && 
            signalLeptons.at(2)->pT() > 20. && 
            nSignalJets >= 5 && signalJets.at(4)->pT() > 30. &&
            // TODO
            // -
            // TODO
            met > 150. &&
            SFOSpair_pTs.at(0) > 150.
            // -
           ) 
          _numSR["SR1B"]++;

        // SR2A
        if (trigger && 
            signalLeptons.at(2)->pT() < 20. && 
            nSignalJets >= 3 && signalJets.at(2)->pT() > 30. && 
            // TODO
            signalJets.at(0)->pT() > 150. &&
            // TODO
            met > 200. &&
            SFOSpair_pTs.at(-1) < 50. 
            // -
           ) 
          _numSR["SR2A"]++;

        // SR2B
        if (trigger && 
            signalLeptons.at(2)->pT() < 60. && 
            nSignalJets >= 3 && signalJets.at(2)->pT() > 30. &&
           // TODO
           // -
           // TODO
           met > 350. && 
           SFOSpair_pTs.at(0) > 150.
           // -
           ) 
          _numSR["SR2B"]++;

      }

      /// Combine the variables of another copy of this analysis (typically on another thread) into this one.
      void combine(const Analysis* other)
      {
        const Analysis_ATLAS_13TeV_2OSLEP_Z_139invfb* specificOther
                = dynamic_cast<const Analysis_ATLAS_13TeV_2OSLEP_Z_139invfb*>(other);

        for (auto& el : _numSR)
        {
          el.second += specificOther->_numSR.at(el.first);
        }

      }

      // This function can be overridden by the derived SR-specific classes
      virtual void collect_results()
      {

        // add_result(SignalRegionData("SR label", n_obs, {s, s_sys}, {b, b_sys}));
        // TODO: Make sure background events and uncertainties are correct
        add_result(SignalRegionData("SR1A", 3., {_numSR["SR1A"], 0.}, {5.4, 0.7}));
        add_result(SignalRegionData("SR1B", 14., {_numSR["SR1B"], 0.}, {12.8, 1.6}));
        add_result(SignalRegionData("SR2A", 3., {_numSR["SR2A"], 0.}, {5.7, 1.7}));
        add_result(SignalRegionData("SR2B", 6., {_numSR["SR2B"], 0.}, {5.4, 0.8}));


      }


    protected:
      void analysis_specific_reset()
      {
        for (auto& el : _numSR) { el.second = 0.;}
        #ifdef CHECK_CUTFLOW
          std::fill(cutFlowVector.begin(), cutFlowVector.end(), 0);
        #endif
      }

    };

    // Factory fn
    DEFINE_ANALYSIS_FACTORY(ATLAS_13TeV_2OSLEP_Z_139invfb)


  }
}
