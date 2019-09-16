//   GAMBIT: Global and Modular BSM Inference Tool
//  *********************************************
///  \file
///
///  Helper functions for converting between
///  different event types.
///
///  *********************************************
///
///  Authors (add name and date if you modify):
///
///  \author Andy Buckley
///  \author Anders Kvellestad
///  \author Pat Scott
///  \author Martin White
///
///  *********************************************

#pragma once

#include "gambit/ColliderBit/colliders/Pythia8/Py8Utils.hpp"
#include "gambit/Elements/shared_types.hpp"//AddedTP Aug19 for debug

#include "HEPUtils/Event.h"
#include "HEPUtils/Particle.h"
#include "HEPUtils/FastJet.h"
#include "MCUtils/PIDCodes.h"

#include <fstream>
#include <algorithm>//Used for sorting the particles so I know CBS and gambit get them in the same order - makes debugging easier.

//#define MET_DEBUG
//#define JET_CHECKING
//#define SINGLE_EVENT_DEBUG

//#define ELECTRON_PARENTAGE_DEBUG
//#define ELECTRON_PARENTAGE_DEBUG_TWO

//#define NJET_DATA_OUTPUT



inline bool compare_particles_by_pz(FJNS::PseudoJet jet1, FJNS::PseudoJet jet2)
{
  return (jet1.pz() > jet2.pz());
}


namespace Gambit
{

  namespace ColliderBit
  {

    #ifndef DODGY_GLOBAL_COUNTING_SHORTCUT
    static int TomeksCounter1{0};
    #define DODGY_GLOBAL_COUNTING_SHORTCUT
    #endif

    /// Convert a hadron-level EventT into an unsmeared HEPUtils::Event
    /// @todo Overlap between jets and prompt containers: need some isolation in MET calculation
    template<typename EventT>
    void convertParticleEvent(const EventT& pevt, HEPUtils::Event& result, double antiktR)
    {
      result.clear();

      //Want to see debug info in more detail.
      //std::cout.precision(6);

      std::vector<FJNS::PseudoJet> bhadrons; //< for input to FastJet b-tagging
      std::vector<HEPUtils::Particle> bpartons, cpartons, tauCandidates;
      HEPUtils::P4 pout; //< Sum of momenta outside acceptance

      // Make a first pass of non-final particles to gather b-hadrons and taus
      for (int i = 0; i < pevt.size(); ++i)
      {
        const auto& p = pevt[i];

        // Find last b-hadrons in b decay chains as the best proxy for b-tagging
        /// @todo Temporarily using quark-based tagging instead -- fix
        if (p.idAbs() == 5) {
          bool isGoodB = true;
          const std::vector<int> bDaughterList = p.daughterList();
          for (size_t daughter = 0; daughter < bDaughterList.size(); daughter++)
          {
            const auto& pDaughter = pevt[bDaughterList[daughter]];
            int daughterID = pDaughter.idAbs();
            if (daughterID == 5) isGoodB = false;
          }
          if (isGoodB)
            bpartons.push_back(HEPUtils::Particle(mk_p4(p.p()), p.id()));
        }

        // Find last c-hadrons in decay chains as the best proxy for c-tagging
        /// @todo Temporarily using quark-based tagging instead -- fix
        if (p.idAbs() == 4) {
          bool isGoodC = true;
          const std::vector<int> cDaughterList = p.daughterList();
          for (size_t daughter = 0; daughter < cDaughterList.size(); daughter++)
          {
            const auto& pDaughter = pevt[cDaughterList[daughter]];
            int daughterID = pDaughter.idAbs();
            if (daughterID == 4) isGoodC = false;
          }
          if (isGoodC)
            cpartons.push_back(HEPUtils::Particle(mk_p4(p.p()), p.id()));
        }

        // Find tau candidates
        if (p.idAbs() == MCUtils::PID::TAU) {
          HEPUtils::P4 tmpMomentum;
          bool isGoodTau=true;
          const std::vector<int> tauDaughterList = p.daughterList();
          for (size_t daughter = 0; daughter < tauDaughterList.size(); daughter++)
          {
            const auto& pDaughter = pevt[tauDaughterList[daughter]];
            int daughterID = pDaughter.idAbs();
            // Veto leptonic taus
            /// @todo What's wrong with having a W daughter? Doesn't that just mark a final tau?
            if (daughterID == MCUtils::PID::ELECTRON || daughterID == MCUtils::PID::MUON ||
                daughterID == MCUtils::PID::WPLUSBOSON || daughterID == MCUtils::PID::TAU)
              isGoodTau = false;
            if (daughterID != MCUtils::PID::TAU) tmpMomentum += mk_p4(pDaughter.p());
          }

          if (isGoodTau) {
            tauCandidates.push_back(HEPUtils::Particle(mk_p4(p.p()), p.id()));
          }
        }
      }

#ifdef SINGLE_EVENT_DEBUG
      std::fstream myfile;
      myfile.open("GAMBIT_SINGLE_EVENT_FILE.csv", std::fstream::app);
      myfile << "Particle ID, Prompt?, Visible?, Px, Py, Pz, Pe, Abs(eta), Parents, Parent Pe, etc.\n";
#endif

#ifdef JET_CHECKING
    std::fstream jetfile;
    if (TomeksCounter1 == 22)
    {  
      jetfile.open("GAMBIT Jetfile.csv", std::fstream::app);
    }
#endif



      // Loop over final state particles for jet inputs and MET
      std::vector<FJNS::PseudoJet> jetparticles;
      for (int i = 0; i < pevt.size(); ++i)
      {
        const auto& p = pevt[i];


        // Only consider final state particles
        if (!p.isFinal()) continue;

        

        // Check there's no partons!!
        if (p.id() == 21 || abs(p.id()) <= 6) {
          std::ostringstream sid;
          bool gotmother = false;
          if (p.mother1() != 0) { gotmother = true; sid << pevt[p.mother1()].id() << " "; }
          if (p.mother2() != 0) { gotmother = true; sid << pevt[p.mother2()].id() << " "; }
          if (gotmother) sid << " -> ";
          sid << p.id();
          ColliderBit_error().forced_throw(LOCAL_INFO, "Found final-state parton " + sid.str() + " in particle-level event converter: "
           "reconfigure your generator to include hadronization, or Gambit to use the partonic event converter.");
        }

        // Add particle outside ATLAS/CMS acceptance to MET
        /// @todo Move out-of-acceptance MET contribution to BuckFast
        if (std::abs(p.eta()) > 5.0)
        {
          pout += mk_p4(p.p());
          #ifdef MET_DEBUG
          std::cout << "Adding particle to met - PID: " << (p.id()) << "; Pz: " << mk_p4(p.p()).pz() << ", m: " << mk_p4(p.p()).m() << std::endl;
          #endif
          continue;
        }

        // Promptness: for leptons and photons we're only interested if they don't come from hadron/tau decays
        const bool prompt = !fromHadron(i, pevt); //&& !fromTau(i, pevt);
        const bool visible = MCUtils::PID::isStrongInteracting(p.id()) || MCUtils::PID::isEMInteracting(p.id());

#ifdef SINGLE_EVENT_DEBUG

        /*std::vector<int>Full_Mother_List;
        int mother_part = i;
        while (pevt[mother_part].mother1() != 0)
        {
           Full_Mother_List.push_back(mother_part);
           mother_part = pevt[mother_part].mother1();
           //if (pevt[mother_part].mother1() != 0)
               //std::cout << "MULTIPLE PARENTS!";
        }
        

        myfile << p.id() << ", " << prompt << ", " << visible << ", " << p.p().px() << ", " << p.p().py() << ", " << p.p().pz() << ", " << p.p().e() << ", " << abs(p.eta()) << ", ";

        for (int parent_part : Full_Mother_List)
        {
            myfile << pevt[parent_part].id() << ", " << pevt[parent_part].e() << ", ";
        }
        myfile << "\n";*/
#endif


        // Add prompt and invisible particles as individual particles
        if (prompt || !visible)
        {
          HEPUtils::Particle* gp = new HEPUtils::Particle(mk_p4(p.p()), p.id());
          gp->set_prompt();
          result.add_particle(gp);
          #ifdef JET_CHECKING
          if (p.id() == 22)
          {
            std::cout << "Prompt Photon: pz = " << p.p().pz() << std::endl;
          }
          #endif
        }

        // All particles other than invisibles and muons are jet constituents
        // Matthias added test to keep non-prompt particles
        if (visible && p.idAbs() != MCUtils::PID::MUON)
        {
           jetparticles.push_back(mk_pseudojet(p.p()));
           FJNS::PseudoJet temp = mk_pseudojet(p.p());
#ifdef JET_CHECKING
        if (TomeksCounter1 == 22)
          {
            jetfile << p.id() << ", " << temp.px() << ", " << temp.py() << ", " << temp.pz() << ", " << temp.e() << ", " << temp.rap() << ",\n";
          }
#endif
        }
        // next case are visible non-prompt muons
        //if (visible && p.idAbs() == MCUtils::PID::MUON && !prompt) jetparticles.push_back(mk_pseudojet(p.p()));
        // next case are non-prompt neutrinos
        //if (!visible && !prompt) jetparticles.push_back(mk_pseudojet(p.p()));
 
      }
#ifdef SINGLE_EVENT_DEBUG
      myfile << "\n\n========\n\n";
      myfile.close();
#endif
#ifdef JET_CHECKING
if (TomeksCounter1 == 22)
{
      jetfile << "\n\n========\n\n";
      jetfile.close();
}
#endif
      /// Jet finding
      /// @todo Choose jet algorithm via detector _settings? Run several algs?
      const FJNS::JetDefinition jet_def(FJNS::antikt_algorithm, antiktR);
      //sort jet particles by p so gambit and CBS have same order.
      std::sort(jetparticles.begin(), jetparticles.end(), compare_particles_by_pz);
      //int TP_TEMP_COUNTER = 0;
      //for (auto particle : jetparticles)
      //{
        //std::cout << "Particle Number: " << TP_TEMP_COUNTER++ << "; Pz: " << particle.pz() << "; Px: " << particle.px() <<std::endl;
      //}

      FJNS::ClusterSequence cseq(jetparticles, jet_def);
      std::vector<FJNS::PseudoJet> pjets = sorted_by_pt(cseq.inclusive_jets(10));

      //std::cout << "\n\nJETCLUSTER_DEBUG_INFO: " << std::endl;
      //std::cout << "AntiktR: " << antiktR << std::endl;
      //std::cout << "pTmin: " << 10 << std::endl;
      //std::cout << "Number of particles passed to algorithm is: " << jetparticles.size() << std::endl;
      //std::cout << "Number of outputjets, unsorted: " << cseq.inclusive_jets(10).size() << std::endl;
      //std::cout << "Number of outputjets, sorted: " << sorted_by_pt(cseq.inclusive_jets(10)).size() << std::endl;


      /// Do jet b-tagging, etc. and add to the Event
      /// @todo Use ghost tagging?
      /// @note We need to _remove_ this b-tag in the detector sim if outside the tracker acceptance!
      for (auto& pj : pjets) {
        HEPUtils::P4 jetMom = HEPUtils::mk_p4(pj);

        /// @todo Replace with HEPUtils::any(bhadrons, [&](const auto& pb){ pj.delta_R(pb) < 0.4 })
        bool isB = false;
        for (HEPUtils::Particle& pb : bpartons) {
          if (jetMom.deltaR_eta(pb.mom()) < 0.4) { ///< @todo Hard-coded radius!!!
            isB = true;
            break;
          }
        }

        bool isC = false;
        for (HEPUtils::Particle& pc : cpartons) {
          if (jetMom.deltaR_eta(pc.mom()) < 0.4) { ///< @todo Hard-coded radius!!!
            isC = true;
            break;
          }
        }

        bool isTau = false;
        for (HEPUtils::Particle& ptau : tauCandidates){
          if (jetMom.deltaR_eta(ptau.mom()) < 0.5){
            isTau = true;
            break;
          }
        }

        // Add to the event (use jet momentum for tau)
        if (isTau) {
          HEPUtils::Particle* gp = new HEPUtils::Particle(HEPUtils::mk_p4(pj), MCUtils::PID::TAU);
          gp->set_prompt();
          result.add_particle(gp);
        }



        result.add_jet(new HEPUtils::Jet(HEPUtils::mk_p4(pj), isB, isC));



      }

      /// Calculate missing momentum
      //
      // From balance of all visible momenta (requires isolation)
      // const std::vector<Particle*> visibles = result.visible_particles();
      // HEPUtils::P4 pvis;
      // for (size_t i = 0; i < visibles.size(); ++i) {
      //   pvis += visibles[i]->mom();
      // }
      // for (size_t i = 0; i < result.jets.size(); ++i) {
      //   pvis += result.jets[i]->mom();
      // }
      // set_missingmom(-pvis);
      //
      // From sum of invisibles, including those out of range
      for (size_t i = 0; i < result.invisible_particles().size(); ++i) {
        pout += result.invisible_particles()[i]->mom();
        #ifdef MET_DEBUG
        auto thingy = result.invisible_particles()[i];
        std::cout << "Adding particle to met - PID: " << thingy->pid() << "; Pz: " << thingy->mom().pz() << ", m: " << thingy->mom().m() << std::endl;
        #endif
      }
      #ifdef MET_DEBUG
      std::cout << "FINAL MISSING MOMENTUM IS: (" << pout.px() <<", "<<pout.py()<<", "<<pout.pz()<<", "<<pout.m()<<")"<<std::endl;
      #endif
      result.set_missingmom(pout);

#ifdef NJET_DATA_OUTPUT
    double sumpT = 0;
    for (int i{0}; i < result.jets().size(); i++)
    {
      sumpT += result.jets()[i]->pT();
    }
    std::fstream njetdatafile;
    njetdatafile.open("GambitNjetFile.csv", std::fstream::app);
    njetdatafile << result.jets().size() << ", " << result.met() << ", " << result.jets()[0]->pT() << ", " << sumpT << "\n";
    njetdatafile.close();
#endif

#ifdef SINGLE_EVENT_DEBUG
++TomeksCounter1;
if (TomeksCounter1 == 3047 || TomeksCounter1 == 3048 || TomeksCounter1 == 3541 || TomeksCounter1 == 3542 || TomeksCounter1 == 3618 || TomeksCounter1 == 3619)
    {
      cout << TomeksCounter1 << ". Gambit Event Information\n";
      cout << "  MET  = " << result.met() << " GeV" << endl;
      cout << "  #e   = " << result.electrons().size() << endl;
      cout << "  #mu  = " << result.muons().size() << endl;
      cout << "  #tau = " << result.taus().size() << endl;
      cout << "  #jet = " << result.jets().size() << endl;
      for (int i{0}; i < result.jets().size(); i++)
      {
        cout << "   pT Jet " << i << ": " << result.jets()[i]->pT() << endl;
      }
      cout << "  #pho  = " << result.photons().size() << endl;
      cout << endl;
    }
#endif

#ifdef JET_CHECKING
      jetfile << "\n\n========\n\n";
      jetfile.close();    
#endif  
    }
  

    /// Convert a partonic (no hadrons) EventT into an unsmeared HEPUtils::Event
    template<typename EventT>
    void convertPartonEvent(const EventT& pevt, HEPUtils::Event& result, double antiktR)
    {
      result.clear();

      std::vector<HEPUtils::Particle> tauCandidates;

      // Make a first pass of non-final particles to gather taus
      for (int i = 0; i < pevt.size(); ++i)
      {
        const auto& p = pevt[i];

        // Find last tau in prompt tau replica chains as a proxy for tau-tagging
        if (p.idAbs() == MCUtils::PID::TAU) {
          std::vector<int> tauDaughterList = p.daughterList();
          HEPUtils::P4 tmpMomentum;
          bool isGoodTau=true;

          for (size_t daughter = 0; daughter < tauDaughterList.size(); daughter++)
          {
            const auto& pDaughter = pevt[tauDaughterList[daughter]];
            int daughterID = pDaughter.idAbs();
            if (daughterID == MCUtils::PID::ELECTRON || daughterID == MCUtils::PID::MUON ||
                daughterID == MCUtils::PID::WPLUSBOSON || daughterID == MCUtils::PID::TAU)
              isGoodTau = false;
            if (daughterID != MCUtils::PID::TAU) tmpMomentum += mk_p4(pDaughter.p());
          }

          if (isGoodTau) {
            tauCandidates.push_back(HEPUtils::Particle(mk_p4(p.p()), p.id()));
          }
        }
      }

      std::vector<FJNS::PseudoJet> jetparticles; //< Pseudojets for input to FastJet
      HEPUtils::P4 pout; //< Sum of momenta outside acceptance

      // Make a single pass over the event to gather final leptons, partons, and photons
      for (int i = 0; i < pevt.size(); ++i)
      {
        const auto& p = pevt[i];

        // We only use "final" partons, i.e. those with no children. So Py8 must have hadronization disabled
        if (!p.isFinal()) continue;

        // Only consider partons within ATLAS/CMS acceptance
        /// @todo We should leave this for the detector sim / analysis to deal with
        if (std::abs(p.eta()) > 5.0) {
          pout += mk_p4(p.p());
          continue;
        }

        // Find electrons/muons/taus/photons to be treated as prompt (+ invisibles)
        /// @todo *Some* photons should be included in jets!!! Ignore for now since no FSR
        /// @todo Lepton dressing
        const bool prompt = isFinalPhoton(i, pevt) || (isFinalLepton(i, pevt)); // && std::abs(p.id()) != MCUtils::PID::TAU);
        const bool visible = MCUtils::PID::isStrongInteracting(p.id()) || MCUtils::PID::isEMInteracting(p.id());
        if (prompt || !visible) {
          HEPUtils::Particle* gp = new HEPUtils::Particle(mk_p4(p.p()), p.id());
          gp->set_prompt();
          result.add_particle(gp);
        }

        // Everything other than invisibles and muons, including taus & partons are jet constituents
        /// @todo Only include hadronic tau fraction?
        // if (visible && (isFinalParton(i, pevt) || isFinalTau(i, pevt))) {
        if (visible && p.idAbs() != MCUtils::PID::MUON) {
          FJNS::PseudoJet pj = mk_pseudojet(p.p());
          //pj.set_user_index(std::abs(p.id()));
          jetparticles.push_back(pj);
          
        }

      }

      /// Jet finding
      /// @todo choose jet algorithm via _settings?
      const FJNS::JetDefinition jet_def(FJNS::antikt_algorithm, antiktR);
      FJNS::ClusterSequence cseq(jetparticles, jet_def);
      std::vector<FJNS::PseudoJet> pjets = sorted_by_pt(cseq.inclusive_jets(10));
      // Add to the event, with b-tagging info"
      for (const FJNS::PseudoJet& pj : pjets) {
        // Do jet b-tagging, etc. by looking for b quark constituents (i.e. user index = |parton ID| = 5)
        /// @note This b-tag is removed in the detector sim if outside the tracker acceptance!
        const bool isB = HEPUtils::any(pj.constituents(),
                 [](const FJNS::PseudoJet& c){ return c.user_index() == MCUtils::PID::BQUARK; });
        const bool isC = HEPUtils::any(pj.constituents(),
                 [](const FJNS::PseudoJet& c){ return c.user_index() == MCUtils::PID::CQUARK; });
        result.add_jet(new HEPUtils::Jet(HEPUtils::mk_p4(pj), isB, isC));

        bool isTau=false;
        for(auto& ptau : tauCandidates){
          HEPUtils::P4 jetMom = HEPUtils::mk_p4(pj);
          if(jetMom.deltaR_eta(ptau.mom()) < 0.5){
            isTau=true;
            break;
          }
        }
        // Add to the event (use jet momentum for tau)
        if (isTau) {
          HEPUtils::Particle* gp = new HEPUtils::Particle(HEPUtils::mk_p4(pj), MCUtils::PID::TAU);
          gp->set_prompt();
          result.add_particle(gp);
        }
      }

      /// Calculate missing momentum
      //
      // From balance of all visible momenta (requires isolation)
      // const std::vector<Particle*> visibles = result.visible_particles();
      // HEPUtils::P4 pvis;
      // for (size_t i = 0; i < visibles.size(); ++i) {
      //   pvis += visibles[i]->mom();
      // }
      // for (size_t i = 0; i < result.jets.size(); ++i) {
      //   pvis += result.jets[i]->mom();
      // }
      // set_missingmom(-pvis);
      //
      // From sum of invisibles, including those out of range
      for (const HEPUtils::Particle* p : result.invisible_particles())
        pout += p->mom();
      result.set_missingmom(pout);
    }

  }

}
