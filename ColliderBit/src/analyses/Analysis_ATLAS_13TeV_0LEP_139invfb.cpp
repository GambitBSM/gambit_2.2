// -*- C++ -*-
#include "gambit/ColliderBit/analyses/Analysis.hpp"
#include "gambit/ColliderBit/analyses/Cutflow.hpp"
#include "gambit/ColliderBit/ATLASEfficiencies.hpp"
#include "Eigen/Eigen"

namespace Gambit {
  namespace ColliderBit {

    using namespace std;
    using namespace HEPUtils;


    /// @brief ATLAS Run 2 0-lepton jet+MET SUSY analysis, with 139/fb of data
    ///
    /// Based on:
    ///   https://cds.cern.ch/record/2686254
    ///   https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/CONFNOTES/ATLAS-CONF-2019-040/
    ///
    class Analysis_ATLAS_13TeV_0LEP_139invfb : public Analysis {
    public:

      // Required detector sim
      static constexpr const char* detector = "ATLAS";

      // Numbers passing cuts
      static const size_t NUMSR = 10;
      double _srnums[NUMSR] = {0., 0., 0., 0., 0., 0., 0., 0., 0., 0.};
      Cutflows _cutflows;
      enum SRNames { SR2J_1600=0, SR2J_2200, SR2J_2800,
                     SR4J_1000, SR4J_2200, SR4J_3400, SR5J_1600,
                     SR6J_1000, SR6J_2200, SR6J_3400 };

      Analysis_ATLAS_13TeV_0LEP_139invfb() {

        set_analysis_name("ATLAS_13TeV_0LEP_139invfb");
        set_luminosity(139.0);

        // Book cut-flows
        const vector<string> cutnames = {"Pre-sel", "Njet", "pT1", "pTx", "|eta_x|", "Dphi(j123,MET)min", "Dphi(j4+,MET)min", "Aplanarity", "MET/sqrt(HT)", "m_eff(incl)"};
        _cutflows.addCutflow("2j-1600", cutnames);
        _cutflows.addCutflow("2j-2200", cutnames);
        _cutflows.addCutflow("2j-2800", cutnames);
        _cutflows.addCutflow("4j-1000", cutnames);
        _cutflows.addCutflow("4j-2200", cutnames);
        _cutflows.addCutflow("4j-3400", cutnames);
        _cutflows.addCutflow("5j-1600", cutnames);
        _cutflows.addCutflow("6j-1000", cutnames);
        _cutflows.addCutflow("6j-2200", cutnames);
        _cutflows.addCutflow("6j-3400", cutnames);

      }

      void run(const Event* event) {

        _cutflows.fillinit();

        // Missing energy
        /// @todo Compute from hard objects instead?
        const P4 pmiss = event->missingmom();
        const double met = event->met();


        // Get baseline jets
        /// @todo Drop b-tag if pT < 50 GeV or |eta| > 2.5?
        vector<const Jet*> baselineJets;
        for (const Jet* jet : event->jets()) {
          if (jet->pT() > 20. && jet->abseta() < 2.8) {
            baselineJets.push_back(jet);
          }
        }


        /// @todo Apply a random 9% loss / 0.91 reweight for jet quality criteria?

        // Get baseline electrons and apply efficiency
        vector<Particle*> baselineElectrons;
        for (Particle* electron : event->electrons()) {
          if (electron->pT() > 7. && electron->abseta() < 2.47)
            baselineElectrons.push_back(electron);
        }
        ATLAS::applyElectronEff(baselineElectrons);

        // Get baseline muons and apply efficiency
        vector<Particle*> baselineMuons;
        for (Particle* muon : event->muons()) {
          if (muon->pT() > 6. && muon->abseta() < 2.7)
            baselineMuons.push_back(muon);
        }
        ATLAS::applyMuonEff(baselineMuons);

        // Remove any |eta| < 2.8 jet within dR = 0.2 of an electron
        vector<const Jet*> signalJets;
        for (const Jet* j : baselineJets)
          if (all_of(baselineElectrons, [&](const Particle* e){ return deltaR_rap(*e, *j) > 0.2; }))
            signalJets.push_back(j);

        // Remove electrons with dR = shrinking cone of surviving |eta| < 2.8 jets
        vector<const Particle*> signalElectrons;
        for (const Particle* e : baselineElectrons)
          if (all_of(signalJets, [&](const Jet* j){ return deltaR_rap(*e, *j) > min(0.4, 0.04+10/e->pT()); }))
            signalElectrons.push_back(e);
        // Apply electron ID selection
        ATLAS::applyLooseIDElectronSelectionR2(signalElectrons);
        /// @todo And tight ID for high purity... used where?

        // Remove muons with dR = 0.4 of surviving |eta| < 2.8 jets
        /// @note Within 0.2, discard the *jet* based on jet track vs. muon criteria... can't be done yet
        vector<const Particle*> signalMuons;
        for (const Particle* m : baselineMuons)
          if (all_of(signalJets, [&](const Jet* j){ return deltaR_rap(*m, *j) > min(0.4, 0.04+10/m->pT()); }))
            signalMuons.push_back(m);
        /// @todo And tight ID for high purity... used where?

        // The subset of jets with pT > 50 GeV is used for several calculations
        vector<const Jet*> signalJets50;
        for (const Jet* j : signalJets)
          if (j->pT() > 50) signalJets50.push_back(j);


        ////////////////////////////////
        // Calculate common variables and cuts

        // Multiplicities
        const size_t nElectrons = signalElectrons.size();
        const size_t nMuons = signalMuons.size();
        const size_t nJets50 = signalJets50.size();
        const size_t nJets = signalJets.size();

        // HT-related quantities (calculated over all >50 GeV jets)
        double sumptj = 0;
        for (const Jet* j : signalJets50) sumptj += j->pT();
        const double HT = sumptj;
        const double sqrtHT = sqrt(HT);
        const double met_sqrtHT = met/sqrtHT;

        // Meff-related quantities (calculated over >50 GeV jets only)
        double sumptj50_incl = 0; // sumptj50_4 = 0, sumptj50_5 = 0, sumptj50_6 = 0;
        for (size_t i = 0; i < signalJets50.size(); ++i) {
          const Jet* j = signalJets50[i];
          // if (i < 4) sumptj50_4 += j->pT();
          // if (i < 5) sumptj50_5 += j->pT();
          // if (i < 6) sumptj50_6 += j->pT();
          sumptj50_incl += j->pT();
        }
        // const double meff_4 = met + sumptj50_4;
        // const double meff_5 = met + sumptj50_5;
        // const double meff_6 = met + sumptj50_6;
        // const double meff_incl = met + sumptj50_incl;
        const double meff = met + sumptj50_incl;
        // const double met_meff_4 = met / meff_4;
        // const double met_meff_5 = met / meff_5;
        // const double met_meff_6 = met / meff_6;

        // Jet |eta|s
        double etamax_2 = 0, etamax_4 = 0, etamax_5 = 0, etamax_6 = 0;
        for (size_t i = 0; i < signalJets50.size(); ++i) {
          const Jet* j = signalJets50[i];
          if (i < 2) etamax_2 = max(etamax_2, j->abseta());
          if (i < 4) etamax_4 = max(etamax_4, j->abseta());
          if (i < 5) etamax_5 = max(etamax_5, j->abseta());
          if (i < 6) etamax_6 = max(etamax_6, j->abseta());
        }

        // Jet--MET dphis
        double dphimin_123 = DBL_MAX, dphimin_more = DBL_MAX;
        for (size_t i = 0; i < min(3lu,signalJets50.size()); ++i)
          dphimin_123 = min(dphimin_123, acos(cos(signalJets50[i]->phi() - pmiss.phi())));
        for (size_t i = 3; i < signalJets50.size(); ++i)
          dphimin_more = min(dphimin_more, acos(cos(signalJets50[i]->phi() - pmiss.phi())));

        // Jet aplanarity (on 50 GeV jets only, cf. paper)
        Eigen::Matrix3d momtensor = Eigen::Matrix3d::Zero();
        double norm = 0;
        for (const Jet* jet : signalJets50) {
          const P4& p4 = jet->mom();
          norm += p4.p2();
          for (size_t i = 0; i < 3; ++i) {
            const double pi = (i == 0) ? p4.px() : (i == 1) ? p4.py() : p4.pz();
            for (size_t j = 0; j < 3; ++j) {
              const double pj = (j == 0) ? p4.px() : (j == 1) ? p4.py() : p4.pz();
              momtensor(i,j) += pi*pj;
            }
          }
        }
        momtensor /= norm;
        const double mineigenvalue = momtensor.eigenvalues().real().minCoeff();
        const double aplanarity = 1.5 * mineigenvalue;


        ////////////////////////////////
        // Fill signal regions and cutflows

        // Preselection
        if (nElectrons + nMuons != 0) return;
        if (nJets50 < 2 || signalJets50[0]->pT() < 200) return;
        if (met < 300) return;
        if (meff < 800) return;
        if (dphimin_123 < 0.4) return;

        _cutflows.fill(0);
        const double w = event->weight();

        // 2 jet regions
        if (nJets >= 2) {
          if (_cutflows["2j-1600"].filltail({
                nJets50 >= 2, signalJets[0]->pT() > 250, signalJets[1]->pT() > 250, etamax_2 < 2.0,
                dphimin_123 > 0.8, dphimin_more > 0.4, true, met_sqrtHT > 16, meff > 1600})) _srnums[SR2J_1600] += w;
          if (_cutflows["2j-2200"].filltail({
                nJets50 >= 2, signalJets[0]->pT() > 600, signalJets[1]->pT() >  50, etamax_2 < 2.8,
                dphimin_123 > 0.4, dphimin_more > 0.2, true, met_sqrtHT > 16, meff > 2200})) _srnums[SR2J_2200] += w;
          if (_cutflows["2j-2800"].filltail({
                nJets50 >= 2, signalJets[0]->pT() > 250, signalJets[1]->pT() > 250, etamax_2 < 1.2,
                dphimin_123 > 0.8, dphimin_more > 0.4, true, met_sqrtHT > 16, meff > 2800})) _srnums[SR2J_2800] += w;
        }

        // 4 jet regions
        if (nJets >= 4) {
          if (_cutflows["4j-1000"].filltail({
                nJets50 >= 4, signalJets.at(0)->pT() > 200, signalJets.at(3)->pT() > 100, etamax_4 < 2.0,
                dphimin_123 > 0.4, dphimin_more > 0.4, true, met_sqrtHT > 16, meff > 1000})) _srnums[SR4J_1000] += w;
          if (_cutflows["4j-2200"].filltail({
                nJets50 >= 4, signalJets[0]->pT() > 200, signalJets[3]->pT() > 100, etamax_4 < 2.0,
                dphimin_123 > 0.4, dphimin_more > 0.4, true, met_sqrtHT > 16, meff > 2200})) _srnums[SR4J_2200] += w;
          if (_cutflows["4j-3400"].filltail({
                nJets50 >= 4, signalJets[0]->pT() > 200, signalJets[3]->pT() > 100, etamax_4 < 2.0,
                dphimin_123 > 0.4, dphimin_more > 0.4, true, met_sqrtHT > 10, meff > 3400})) _srnums[SR4J_3400] += w;
        }

        // 5 jet region
        if (nJets >= 5) {
          if (_cutflows["5j-1600"].filltail({
                nJets50 >= 5, signalJets[0]->pT() > 600, signalJets[4]->pT() > 50, etamax_5 < 2.8,
                dphimin_123 > 0.4, dphimin_more > 0.2, true, met_sqrtHT > 16, meff > 1600})) _srnums[SR5J_1600] += w;
        }

        // 6 jet regions
        if (nJets >= 6) {
          if (_cutflows["6j-1000"].filltail({
                nJets50 >= 6, signalJets[0]->pT() > 200, signalJets[5]->pT() > 75, etamax_6 < 2.0,
                dphimin_123 > 0.4, dphimin_more > 0.2, aplanarity > 0.08, met_sqrtHT > 16, meff > 1000})) _srnums[SR6J_1000] += w;
          if (_cutflows["6j-2200"].filltail({
                nJets50 >= 6, signalJets[0]->pT() > 200, signalJets[5]->pT() > 75, etamax_6 < 2.0,
                dphimin_123 > 0.4, dphimin_more > 0.2, aplanarity > 0.08, met_sqrtHT > 16, meff > 2200})) _srnums[SR6J_2200] += w;
          if (_cutflows["6j-3400"].filltail({
                nJets50 >= 6, signalJets[0]->pT() > 200, signalJets[5]->pT() > 75, etamax_6 < 2.0,
                dphimin_123 > 0.4, dphimin_more > 0.2, aplanarity > 0.08, met_sqrtHT > 10, meff > 3400})) _srnums[SR6J_3400] += w;
        }


        // if (dphimin_123 > 0.8 && dphimin_more > 0.4) {
        //   if (signalJets[1]->pT() > 250 && etamax_2 < 0.8) { //< implicit pT[0] cut
        //     if (met_sqrtHT > 14 && meff_incl > 1200) num_2j_1200 += 1;
        //   }
        //   if (signalJets[1]->pT() > 300 && etamax_2 < 1.2) { //< implicit pT[0] cut
        //     if (met_sqrtHT > 18 && meff_incl > 1600) num_2j_1600 += 1;
        //   }
        //   if (signalJets[1]->pT() > 350 && etamax_2 < 1.2) { //< implicit pT[0] cut
        //     if (met_sqrtHT > 18 && meff_incl > 2000) num_2j_2000 += 1;
        //     if (met_sqrtHT > 18 && meff_incl > 2400) num_2j_2400 += 1;
        //     if (met_sqrtHT > 18 && meff_incl > 2800) num_2j_2800 += 1;
        //   }
        //   if (signalJets[1]->pT() > 350) { //< implicit pT[0] cut
        //     if (met_sqrtHT > 18 && meff_incl > 3600) num_2j_3600 += 1;
        //   }
        // }

        // if (dphimin_123 > 0.4 && dphimin_more > 0.2) {
        //   if(signalJets[0]->pT() > 600 && signalJets[1]->pT() > 50){
        //     if (met_sqrtHT > 26 && meff_incl > 2100) num_2j_2100 += 1;
        //   }
        // }

        // // 4 jet regions (note implicit pT[1,2] cuts)
        // if (nJets50 >= 4 && dphimin_123 > 0.4 && dphimin_more > 0.4 && signalJets[0]->pT() > 200 && aplanarity > 0.04) {
        //   if (signalJets[3]->pT() > 100 && etamax_4 < 1.2 && met_meff_4 > 0.3 && meff_incl > 1000) num_4j_1000 += 1;
        //   if (signalJets[3]->pT() > 100 && etamax_4 < 2.0 && met_meff_4 > 0.25 && meff_incl > 1400) num_4j_1400 += 1;
        //   if (signalJets[3]->pT() > 100 && etamax_4 < 2.0 && met_meff_4 > 0.25 && meff_incl > 1800) num_4j_1800 += 1;
        //   if (signalJets[3]->pT() > 100 && etamax_4 < 2.0 && met_meff_4 > 0.25 && meff_incl > 2200) num_4j_2200 += 1;
        //   if (signalJets[3]->pT() > 150 && etamax_4 < 2.0 && met_meff_4 > 0.20 && meff_incl > 2600) num_4j_2600 += 1;
        //   if (signalJets[3]->pT() > 150 && etamax_4 < 2.0 && met_meff_4 > 0.20 && meff_incl > 3000) num_4j_3000 += 1;
        // }

        // // 5 jet regions (note implicit pT[1,2,3] cuts)
        // if (nJets50 >= 5){

        //   if(signalJets[0]->pT() > 700. && signalJets[4]->pT() > 50. && dphimin_123 > 0.4 && dphimin_more > 0.2 && met_meff_5 > 0.3 &&  meff_incl > 1700) num_5j_1700 += 1;
        //   if(signalJets[0]->pT() > 200. && signalJets[4]->pT() > 50. && dphimin_123 > 0.4 && dphimin_more > 0.2 && met_meff_5 > 0.15 &&  aplanarity > 0.08 && meff_incl > 1600) num_5j_1600 += 1;
        //   if(signalJets[0]->pT() > 200. && signalJets[4]->pT() > 50. && dphimin_123 > 0.4 && dphimin_more > 0.4 && met_sqrtHT > 15 && meff_incl > 2000) num_5j_2000 += 1;
        //   if(signalJets[0]->pT() > 200. && signalJets[4]->pT() > 50. && dphimin_123 > 0.8 && dphimin_more > 0.4 && met_sqrtHT > 18 && meff_incl > 2600) num_5j_2600 += 1;

        // }

        // // 6 jet regions (note implicit pT[1,2,3,4] cuts)
        // if (nJets50 >= 6 && dphimin_123 > 0.4 && dphimin_more > 0.2 && signalJets[0]->pT() > 200) {
        //   if (signalJets[5]->pT() >  50 && etamax_6 < 2.0 && met_meff_6 > 0.25 && meff_incl > 1200) num_6j_1200 += 1;
        //   if (signalJets[5]->pT() > 100 && etamax_6 < 2.0 && met_meff_6 > 0.2 && aplanarity > 0.04 && meff_incl > 1800) num_6j_1800 += 1;
        //   if (signalJets[5]->pT() > 100 &&                   met_meff_6 > 0.2 && aplanarity > 0.08 && meff_incl > 2200) num_6j_2200 += 1;
        //   if (signalJets[5]->pT() > 100 &&                   met_meff_6 > 0.15 && aplanarity > 0.08 && meff_incl > 2600) num_6j_2600 += 1;
        // }

        // // Cutflows
        // const vector<string> cuts23j = {"Pre-sel+MET+pT1+meff", "Njet", "Dphi_min(j123,MET)", "Dphi_min(j4+,MET)", "pT2", "eta_j12", "MET/sqrtHT", "m_eff(incl)"};

        // if (nJets >= 2) _flows["2j-1200"].filltail({meff_incl > 800 && signalJets[0]->pT() > 200. , dphimin_123 > 0.8, dphimin_more > 0.4, signalJets[1]->pT() > 250, etamax_2 < 0.8, met_sqrtHT > 14, meff_incl > 1200});
        // if (nJets >= 2) _flows["2j-1600"].filltail({meff_incl > 800 && signalJets[0]->pT() > 200., dphimin_123 > 0.8, dphimin_more > 0.4, signalJets[1]->pT() > 300, etamax_2 < 1.2, met_sqrtHT > 18, meff_incl > 1600});
        // if (nJets >= 2) _flows["2j-2000"].filltail({meff_incl > 800 && signalJets[0]->pT() > 200., dphimin_123 > 0.8, dphimin_more > 0.4, signalJets[1]->pT() > 350, etamax_2 < 1.2, met_sqrtHT > 18, meff_incl > 2000});
        // if (nJets >= 2) _flows["2j-2100"].filltail({meff_incl > 800 && signalJets[0]->pT() > 200., dphimin_123 > 0.4, dphimin_more > 0.2, signalJets[0]->pT() > 600, true, met_sqrtHT > 26, meff_incl > 2100});
        // if (nJets >= 2) _flows["2j-2400"].filltail({meff_incl > 800 && signalJets[0]->pT() > 200., dphimin_123 > 0.8, dphimin_more > 0.4, signalJets[1]->pT() > 350, etamax_2 < 1.2, met_sqrtHT > 18, meff_incl > 2400});
        // if (nJets >= 2) _flows["2j-2800"].filltail({meff_incl > 800 && signalJets[0]->pT() > 200., dphimin_123 > 0.8, dphimin_more > 0.4, signalJets[1]->pT() > 350, etamax_2 < 1.2, met_sqrtHT > 18, meff_incl > 2800});
        // if (nJets >= 2) _flows["2j-3600"].filltail({meff_incl > 800 && signalJets[0]->pT() > 200., dphimin_123 > 0.8, dphimin_more > 0.4, signalJets[1]->pT() > 350, true, met_sqrtHT > 18, meff_incl > 3600});

        // if (nJets >= 3) _flows["3j-1300"].filltail({meff_incl > 800 && signalJets[0]->pT() > 200. , dphimin_123 > 0.4, dphimin_more > 0.2, signalJets[0]->pT() > 700, true, met_sqrtHT > 18, meff_incl > 1300});

        // //const vector<string> cuts456j = {"Pre-sel+MET+pT1+meff", "Njet", "Dphi_min(j123,MET)", "Dphi_min(j4+,MET)", "pT4", "eta_j1234", "Aplanarity", "MET/m_eff(Nj)", "m_eff(incl)"};
        // if (nJets >= 4) _flows["4j-1000"].filltail({meff_incl > 800 && signalJets[0]->pT() > 200. , nJets>=4, dphimin_123 > 0.4, dphimin_more > 0.4, signalJets[3]->pT() > 100, etamax_4 < 1.2, aplanarity > 0.04, met_meff_4 > 0.3, meff_incl > 1000});
        // if (nJets >= 4) _flows["4j-1400"].filltail({meff_incl > 800 && signalJets[0]->pT() > 200. , nJets>=4, dphimin_123 > 0.4, dphimin_more > 0.4, signalJets[3]->pT() > 100, etamax_4 < 2.0, aplanarity > 0.04, met_meff_4 > 0.25, meff_incl > 1400});
        // if (nJets >= 4) _flows["4j-1800"].filltail({meff_incl > 800 && signalJets[0]->pT() > 200. , nJets>=4, dphimin_123 > 0.4, dphimin_more > 0.4, signalJets[3]->pT() > 100, etamax_4 < 2.0, aplanarity > 0.04, met_meff_4 > 0.25, meff_incl > 1800});
        // if (nJets >= 4) _flows["4j-2200"].filltail({meff_incl > 800 && signalJets[0]->pT() > 200. , nJets>=4, dphimin_123 > 0.4, dphimin_more > 0.4, signalJets[3]->pT() > 100, etamax_4 < 2.0, aplanarity > 0.04, met_meff_4 > 0.25, meff_incl > 2200});
        // if (nJets >= 4) _flows["4j-2600"].filltail({meff_incl > 800 && signalJets[0]->pT() > 200. , nJets>=4, dphimin_123 > 0.4, dphimin_more > 0.4, signalJets[3]->pT() > 150, etamax_4 < 2.0, aplanarity > 0.04, met_meff_4 > 0.2, meff_incl > 2600});
        // if (nJets >= 4) _flows["4j-3000"].filltail({meff_incl > 800 && signalJets[0]->pT() > 200. , nJets>=4, dphimin_123 > 0.4, dphimin_more > 0.4, signalJets[3]->pT() > 150, etamax_4 < 2.0, aplanarity > 0.04, met_meff_4 > 0.2, meff_incl > 3000});

        // if (nJets >= 5) _flows["5j-1600"].filltail({meff_incl > 800 && signalJets[0]->pT() > 200., nJets>=5, dphimin_123 > 0.4, dphimin_more > 0.2, true, true, aplanarity > 0.08, met_meff_5 > 0.15, meff_incl > 1600});
        // if (nJets >= 5) _flows["5j-1700"].filltail({meff_incl > 800 && signalJets[0]->pT() > 200., nJets>=5, dphimin_123 > 0.4, dphimin_more > 0.2, signalJets[0]->pT() > 700., true, true, met_meff_5 > 0.3, meff_incl > 1700});
        // if (nJets >= 6) _flows["6j-1200"].filltail({meff_incl > 800 && signalJets[0]->pT() > 200., nJets>=6, dphimin_123 > 0.4, dphimin_more > 0.2, true, etamax_6 < 2.0, true, met_meff_6 > 0.25, meff_incl > 1200});
        // if (nJets >= 6) _flows["6j-1800"].filltail({meff_incl > 800 && signalJets[0]->pT() > 200., nJets>=6, dphimin_123 > 0.4, dphimin_more > 0.2, signalJets[5]->pT() > 100, etamax_6 < 2.0, aplanarity > 0.04, met_meff_6 > 0.2, meff_incl > 1800});
        // if (nJets >= 6) _flows["6j-2200"].filltail({meff_incl > 800 && signalJets[0]->pT() > 200., nJets>=6, dphimin_123 > 0.4, dphimin_more > 0.2, signalJets[5]->pT() > 100, true, aplanarity > 0.08, met_meff_6 > 0.2, meff_incl > 2200});
        // if (nJets >= 6) _flows["6j-2600"].filltail({meff_incl > 800 && signalJets[0]->pT() > 200., nJets>=6, dphimin_123 > 0.4, dphimin_more > 0.2, signalJets[5]->pT() > 100, true, aplanarity > 0.08, met_meff_6 > 0.15, meff_incl > 2600});

      }


      /// Combine the variables of another copy of this analysis (typically on another thread) into this one.
      void combine(const Analysis* other)
      {
        const Analysis_ATLAS_13TeV_0LEP_139invfb* specificOther = dynamic_cast<const Analysis_ATLAS_13TeV_0LEP_139invfb*>(other);
        for (size_t i = 0; i < NUMSR; ++i) _srnums[i] += specificOther->_srnums[i];
      }


      /// Register results objects with the results for each SR; obs & bkg numbers from the CONF note
      void collect_results() {
        add_result(SignalRegionData("SR-2j-1600", 2111, {_srnums[SR2J_1600], 0.}, {2190., 130.}));
        add_result(SignalRegionData("SR-2j-2200",  971, {_srnums[SR2J_2200], 0.}, { 980.,  50.}));
        add_result(SignalRegionData("SR-2j-2800",   78, {_srnums[SR2J_2800], 0.}, {  87.,   8.}));
        add_result(SignalRegionData("SR-4j-1000",  535, {_srnums[SR4J_1000], 0.}, { 536.,  31.}));
        add_result(SignalRegionData("SR-4j-2200",   60, {_srnums[SR4J_2200], 0.}, {  60.,   5.}));
        add_result(SignalRegionData("SR-4j-3400",    4, {_srnums[SR4J_3400], 0.}, {  5.7,  1.0}));
        add_result(SignalRegionData("SR-5j-1600",  320, {_srnums[SR5J_1600], 0.}, { 319.,  19.}));
        add_result(SignalRegionData("SR-6j-1000",   25, {_srnums[SR6J_1000], 0.}, {  21.,  2.9}));
        add_result(SignalRegionData("SR-6j-2200",    5, {_srnums[SR6J_2200], 0.}, {  4.6,  1.0}));
        add_result(SignalRegionData("SR-6j-3400",    0, {_srnums[SR6J_3400], 0.}, {  0.8,  0.4}));
        // const double sf = 139*crossSection()/femtobarn/sumOfWeights();
        // _cutflows.scale(sf);
        // cout << "CUTFLOWS:\n\n" << _cutflows << endl;
      }


    protected:
      void analysis_specific_reset() {
        for (size_t i = 0; i < NUMSR; ++i) _srnums[i] = 0;
      }



    };


    // Factory fn
    DEFINE_ANALYSIS_FACTORY(ATLAS_13TeV_0LEP_139invfb)


  }
}
