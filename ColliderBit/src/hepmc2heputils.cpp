// -*- C++ -*-
//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///
///  hepmc2heputils: a HepMC -> HEPUtils::Event converter
///
///  *********************************************
///
///  Authors (add name and date if you modify): 
///
///  \author Andy Buckley
///  (andy.buckley@cern.ch)
///  \date June 2019
///
///  \author Anders Kvellestad
///  (anders.kvellestad@imperial.ac.uk)
///  \date June 2019
///
///  \author Tomasz Procter
///  (tsp116@ic.ac.uk)
///  \date July-Aug 2019
///  *********************************************

#include "gambit/cmake/cmake_variables.hpp"

#ifndef EXCLUDE_HEPMC

using namespace std;

#include <iostream>

#include "gambit/ColliderBit/hepmc2heputils.hpp"

#include "HepMC3/GenEvent.h"
#include "HepMC3/GenParticle.h"
#include "HEPUtils/FastJet.h"

#include "MCUtils/PIDUtils.h"//Added Tomek Procter July 2019
#include <fstream>
#include <algorithm>//Used for sorting the particles so I know CBS and gambit get them in the same order.

//#define COLLIDERBIT_DEBUG

using namespace HEPUtils;
using namespace FJNS;

//Tomek Procter Aug 19. Note the MCUtils isParton function only checks
// for quarks/gluons (while in pythia, the function used in Gambit inlcudes
//diquarks too), so I'm manually defining this function using the isParton
//is diquark options in MCUtils
inline bool isParton(int pid)
{
   return (MCUtils::PID::isParton(pid) || MCUtils::PID::isDiquark(pid));
}

//TP Aug 19 - Adapted from Pythia function of same name.
inline bool fromHadron(HepMC3::ConstGenParticlePtr gp, int counter = 0)
{
      if (MCUtils::PID::isHadron(abs(gp->pid()))) return true;
      if (isParton(abs(gp->pid()))) return false; // stop the walking at the end of the hadron level
      auto parent_vector = (gp->parents());
      if (parent_vector.size() == 0) return false;
      for (const HepMC3::ConstGenParticlePtr parent : parent_vector)
      {
        if (counter >= 100) return false;//Probably don't need the check to make sure we don't get stuck in an infinite loop.
        if (fromHadron(parent, ++counter)) return true;
      }
      return false;
}


//Tomek Procter Aug 19 - Using the function in MCUtils/HepMCVectors.h, seems
//to cause compile errors. Have implemented this for now, ideally it'd be replaced.
inline double get_eta(HepMC3::ConstGenParticlePtr gp)
{  
    const HepMC3::FourVector& hp4 = gp->momentum();
    double magnitude = sqrt(hp4.px()*hp4.px() + hp4.py()*hp4.py() + hp4.pz()*hp4.pz());
    return atanh(hp4.pz()/magnitude);
}

/// Extract a HepMC event as a HEPUtils::Event
void get_HEPUtils_event(const HepMC3::GenEvent& ge, HEPUtils::Event& evt)
{
  P4 vmet;
  vector<PseudoJet> jetparticles;
  evt.set_weight(ge.weight());

  //Tomek Procter Aug 2019: trying to add B, C, tau checking as in Py8EventConversions.hpp
  std::vector<HEPUtils::Particle> bpartons, cpartons, tauCandidates;

  // Loop over all particles in the event
  for (HepMC3::ConstGenParticlePtr gp : ge.particles())
  {
    // Get status and PID code
    const int st = gp->status();
    const int pid = gp->pid();
    const int apid = fabs(pid);

    // Get 4-momentum
    const HepMC3::FourVector& hp4 = gp->momentum();

    //MODIFIED SO IT NOW USES mkXZE not mkXYZM !!! Tomek Procter Aug 2019
    const P4 p4 = P4::mkXYZE(hp4.px(), hp4.py(), hp4.pz(), hp4.e());

    //---------------------------------------------------------------------------------------------
    //Now we look for b's, c's and tau's a la gambit - Tomek Procter Aug 2019 (Code copied from
    //Py8EventConversions and edited for HepMC where needs be). Would probably ideally be integrated
    //with the gambit code somehow to keep up with changes. The 'fix' comments are copied originally
    //from gambit.

    // Find last b-hadrons in b decay chains as the best proxy for b-tagging
    /// @todo Temporarily using quark-based tagging instead -- fix
    if (apid == 5)
    {
      bool isGoodB = true;
      auto bchildren_vector = gp->children();
      for (const HepMC3::ConstGenParticlePtr bchild : bchildren_vector)
      {
        int childID = abs(bchild->pid());
        if (childID == 5) isGoodB = false;
      }
      if (isGoodB)
      {
        Particle* p = new Particle(p4, pid);
        bpartons.push_back(p);
      }
    }

    // Find last c-hadrons in decay chains as the best proxy for c-tagging
    /// @todo Temporarily using quark-based tagging instead -- fix
    if (apid == 4)
    {
      bool isGoodC = true;
      auto cchildren_vector = gp->children();
      for (const HepMC3::ConstGenParticlePtr cchild : cchildren_vector)
      {
        int daughterID = abs(cchild->pid());
        if (daughterID == 4) isGoodC = false;
      }
      if (isGoodC)
      {
        Particle* p = new Particle(p4, pid);
        cpartons.push_back(p);
      }
    }

    // Find tau candidates
    //Tomek Procter Aug 2019 - I've got rid of the tmpmomentum variable from gambit as it was declared
    //and increased, but nothing is actually done with it 
    if (apid == MCUtils::PID::TAU)
    {
      bool isGoodTau=true;
      auto tauChildList = gp->children();;
      for (const HepMC3::ConstGenParticlePtr tauChild : tauChildList)
      {
        int daughterID = abs(tauChild->pid());
        #ifdef SINGLE_EVENT_DEBUG
        if (TomeksCounter2 >= 499)
        {
          std::cout << "PID: " << daughterID << ", " << "ZMomentum: " << tauChild->momentum().pz() <<std::endl;
        }
        #endif
        // Veto leptonic taus
        /// @todo What's wrong with having a W daughter? Doesn't that just mark a final tau?
        if (daughterID == MCUtils::PID::ELECTRON || daughterID == MCUtils::PID::MUON ||
            daughterID == MCUtils::PID::WPLUSBOSON || daughterID == MCUtils::PID::TAU)
          {
            isGoodTau = false;
          }
      }
      if (isGoodTau)
      {
        Particle* p = new Particle(p4, pid);
        tauCandidates.push_back(p);
        //std::cout << "Is good Tau!" <<std::endl;
      }
    }
    //---------------------------------------------------------------------------------------------
    //END OF B,C,TAU code from Gambit.

    // Use physical particles only
    if (st != 1) continue;//(Note only status code 1 is final - NOT status code 2)

//Here Gambit double checks there are no partons, but not including this hasn't caused any discrepancies.

    //Eta Cut, like in Gambit, for particles outside of detector acceptance.
    if (abs(get_eta(gp)) > 5.0)
    {
      vmet += p4;//Add the momentum of particles outside of acceptance to the MET.
      continue; 
    }
       


// Promptness: for leptons and photons we're only interested if they don't come from hadron/tau decays
    const bool prompt = !fromHadron(gp); //&& !fromTau(i, pevt);<-This was commented out in the original gambit code too.
        

    // Store interacting prompt particles
    /// @todo Dress leptons?
    //Edited Tomek Procter Aug 2019 - Only stores prompt leptons.
    //n.b. the ordering in this section is slightly different to gambit
    //but the logic is the same.
    if ((apid == 22 || apid == 11 || apid == 13 || apid == 15) && prompt)
    {
      Particle* p = new Particle(p4, pid); // the event will take ownership of this pointer
      p->set_prompt(true);
      evt.add_particle(p);
    }
    // Aggregate missing ET
    else if (apid == 12 || apid == 14 || apid == 16 || apid == 1000022 || apid == 1000039)
    {
      Particle* p = new Particle(p4, pid); // the event will take ownership of this pointer
      p->set_prompt(true);
      evt.add_particle(p);
      vmet += p4;
    }
    // Add particles to jets. Note this does now include prompt leptons/photons - 
    if ((apid != MCUtils::PID::MUON) && !(apid == 12 || apid == 14 || apid == 16 || apid == 1000022 || apid == 1000039)) //Tomek Procter August 2019 - Muons don't go into jets, right? So what happens to non-prompt muons?
    {//Now added exclusion of invisibles as n longer else if as we want prompt gammas/leptons.
      PseudoJet pj = PseudoJet(hp4.px(), hp4.py(), hp4.pz(), hp4.e());
      jetparticles.push_back(pj);   
    }
  }

  // MET
  vmet.setPz(0);
  vmet.setM(0);
  evt.set_missingmom(vmet);
  

  // Jets
  vector<PseudoJet> jets = get_jets(jetparticles, 0.4, 10.0); //Cut off momentum edited to match gambit at 10.0
  

  for (const PseudoJet& pj : jets)
  {
    bool hasC = false, hasB = false, hastau = false;
    HEPUtils::P4 jetMom = HEPUtils::mk_p4(pj);
    //Tomek Procter Aug 2019 - replacing the block below with the logic from Gambit.
    //All @todo's are original from Gambit.
    //////////////////////////////////////////////////////////////////////////////////
    /// @todo Bug in HEPUtils::get_jets means that constituent info is lost for now...
    // for (const PseudoJet& c : pj.constituents()) {
    //   if (c.user_index() == 4) hasC = true;
    //   if (c.user_index() == 5) hasB = true;
    // }
    //////////////////////////////////////////////////////////////////////////////////
    for (HEPUtils::Particle& pb : bpartons)
    {
    if (jetMom.deltaR_eta(pb.mom()) < 0.4)///< @todo Hard-coded radius!!!
      {
        hasB = true;
        break;
      }
    }
    for (HEPUtils::Particle& pc : cpartons)
    {
      if (jetMom.deltaR_eta(pc.mom()) < 0.4) { ///< @todo Hard-coded radius!!!
        hasC = true;
        break;
      }
    }
    for (HEPUtils::Particle& ptau : tauCandidates)
    {
      if (jetMom.deltaR_eta(ptau.mom()) < 0.5)
      {
        hastau = true;
        break;
      }
    }
    if (hastau)
    {
      HEPUtils::Particle* tauparticle = new HEPUtils::Particle(HEPUtils::mk_p4(pj), MCUtils::PID::TAU);
      tauparticle->set_prompt();
      evt.add_particle(tauparticle);
    }
    evt.add_jet(new Jet(mk_p4(pj), hasB, hasC));
  }


#ifdef COLLIDERBIT_DEBUG
  // Print event summary
  cout << "\nCBS Event Information\n";
  cout << "  MET  = " << evt.met() << " GeV" << endl;
  cout << "  #e   = " << evt.electrons().size() << endl;
  cout << "  #mu  = " << evt.muons().size() << endl;
  cout << "  #tau = " << evt.taus().size() << endl;
  cout << "  #jet = " << evt.jets().size() << endl;
  for (int i{0}; i < evt.jets().size(); i++)
  {
      cout << "   pT Jet " << i << ": " << evt.jets()[i]->pT() << endl;
  }
  cout << "  #pho  = " << evt.photons().size() << endl;
  cout << endl;
#endif

}

#endif
