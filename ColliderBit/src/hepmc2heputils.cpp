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

int TomeksCounter2{0};

//#define NJET_DATA_OUTPUT

//#define SINGLE_EVENT_DEBUG
//#define ELECTRON_PARENTAGE_DEBUG
//#define JET_CHECKING
//#define MET_DEBUG

//Tomek Procter Aug 19. Note the MCUtils isParton function only checks
// for quarks/gluons (while in pythia, the function used in Gambit inlcudes
//diquarks too), so I'm manually defining this function using the isParton
//is diquark options in MCUtils
inline bool isParton(int pid)
{
   return (MCUtils::PID::isParton(pid) || MCUtils::PID::isDiquark(pid));
}

//#define NJET_DATA_INPUT
#ifdef NJET_DATA_INPUT
std::vector<int> get_GAMBIT_njets_from_CSV()
{
  std::vector<int> to_return;
  to_return.reserve(10000);
  std::fstream myfile;
  myfile.open("GambitNjetFile.csv");
  std::string line;
  while(std::getline(myfile, line))
  {
    to_return.push_back(stoi(line));
    line.clear();
  }
  return to_return;
}
std::vector<int> gambit_njets = get_GAMBIT_njets_from_CSV();
#endif

//TP Aug 19 - Adapted from Pythia function of same name.
inline bool fromHadron(HepMC3::ConstGenParticlePtr gp, bool print_progress = false, int counter = 0)
{
      if (MCUtils::PID::isHadron(abs(gp->pid()))) return true;
      if (isParton(abs(gp->pid()))) return false; // stop the walking at the end of the hadron level
      auto parent_vector = (gp->parents());
      if (parent_vector.size() == 0) return false;
      for (const HepMC3::ConstGenParticlePtr parent : parent_vector)
      {
        if (counter >= 100) return false;
        if (fromHadron(parent, print_progress, ++counter)) return true;
      }
      return false;
}

inline bool compare_particles_by_pz(PseudoJet jet1, PseudoJet jet2)
{
  return (jet1.pz() > jet2.pz());
}


//Tomek Procter Aug 19 - Basically this function's in MCUtils/HepMCVectors.h,
//Except that header file seems a bit broken - including it breaks everything!
//Some of the paths a re funny, think it may have been left befind a bit.
//So here's a manual eta function ideally I'd replace it soon.
inline double get_eta(HepMC3::ConstGenParticlePtr gp)
{  
    const HepMC3::FourVector& hp4 = gp->momentum();
    double magnitude = sqrt(hp4.px()*hp4.px() + hp4.py()*hp4.py() + hp4.pz()*hp4.pz());
    return atanh(hp4.pz()/magnitude);
}

int print_history(HepMC3::ConstGenParticlePtr gp)
{
     std::cout << "Generation 1: ";
     auto parent_vector = (gp->parents());
     if (parent_vector.size() == 0) return 0;
     for (auto parent : parent_vector)
     {
        std::cout << "(" << parent->pid() << ", " << parent->momentum().e() << "), ";
     }
     std::cout << "\nGeneration 2: ";
     for (auto parent : parent_vector)
     {
        auto grandparent_vector = (parent->parents());
        if (grandparent_vector.size() == 0) continue;
        for (auto grandparent : grandparent_vector)
        {
           std::cout << "(" << grandparent->pid() << ", " << grandparent->momentum().e() << "), ";
        }
     }
     std::cout << "\n\n";
     return 0;
}


/// Extract a HepMC event as a HEPUtils::Event
void get_HEPUtils_event(const HepMC3::GenEvent& ge, HEPUtils::Event& evt)
{
  //std::cout << "\n -- -- Begin Event -- -- \n";

  P4 vmet;
  vector<PseudoJet> jetparticles;
  evt.set_weight(ge.weight());

  //Tomek Procter Aug 2019: trying to add B, C, tau checking as in Py8EventConversions.hpp
  std::vector<HEPUtils::Particle> bpartons, cpartons, tauCandidates;


#ifdef SINGLE_EVENT_DEBUG
   std::fstream myfile;
   myfile.open("CBS_SINGLE_EVENT_FILE.csv", std::fstream::app);
   myfile << "Particle ID, Prompt?, Visible?, Px, Py, Pz, Pe, Abs(Eta), Parents, Parent Pe, etc.\n";
#endif

#ifdef JET_CHECKING
std::fstream jetfile;
if (TomeksCounter2 == 22)
{
  jetfile.open("CBS_JetFile.csv", std::fstream::app);
}
#endif

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
    //Now we look for b's, c's and tau's a la gambit - TP Aug 2019 (Code copied from
    //Py8EventConversions and edited where needs be)

    // Find last b-hadrons in b decay chains as the best proxy for b-tagging
        /// @todo Temporarily using quark-based tagging instead -- fix
        if (apid == 5)
        {
          bool isGoodB = true;
          //const std::vector<int> bDaughterList = p.daughterList();
          auto bchildren_vector = gp->children();
          //for (size_t daughter = 0; daughter < bDaughterList.size(); daughter++)
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
          //const std::vector<int> cDaughterList = p.daughterList();
          auto cchildren_vector = gp->children();
          //for (size_t daughter = 0; daughter < cDaughterList.size(); daughter++)
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
        if (apid == MCUtils::PID::TAU)
        {
          P4 tmpMomentum;
          bool isGoodTau=true;
          //const std::vector<int> tauDaughterList = p.daughterList();
          auto tauChildList = gp->children();
          //std::cout << "The Tau Child List is: " <<std::endl;
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
            //if (daughterID != MCUtils::PID::TAU) tmpMomentum += mk_p4(tauChild->momentum());
            //Tomek Procter Aug 2019 - what does tmpMomentum actually do?? Its not used is it? I've commented
            //it out as it seems a waste of time to try and convert it into the right form of 4 momentum.
          }

          if (isGoodTau) {
            Particle* p = new Particle(p4, pid);
            tauCandidates.push_back(p);
            //std::cout << "Is good Tau!" <<std::endl;
          }
        }
      




    //---------------------------------------------------------------------------------------------

    // Use physical particles only
    if (st != 1) continue;

    
    

    //Lets test out adding an eta cut, like they have in gambit:
    if (abs(get_eta(gp)) > 5.0)
    {
    vmet += p4;
    #ifdef MET_DEBUG
    std::cout << "Adding particle to met - PID: " << pid << "; Pz: " << p4.pz() << ", m: " <<p4.m() <<std::endl;
    #endif
    continue; 
    }
       


/////////////////////////////////



// Promptness: for leptons and photons we're only interested if they don't come from hadron/tau decays: Added Tomek Procter Aug 19
        const bool prompt = !fromHadron(gp, (apid == 11)); //&& !fromTau(i, pevt);
        //const bool visible = MCUtils::PID::isStrongInteracting(apid) || MCUtils::PID::isEMInteracting(apid); //Not needed by CBS (yet)

 
#ifdef SINGLE_EVENT_DEBUG
//Find the full parent vector:
/*
    vector<HepMC3::ConstGenParticlePtr> Full_Parent_List;
    if (gp->parent_event() != NULL)
    {
        Full_Parent_List.push_back(gp->parents()[0]);
        int counter{0};
        while (Full_Parent_List.back()->parent_event() != NULL && counter++ < 250 && Full_Parent_List.back()->parents().size() > 0)
        {
            Full_Parent_List.push_back(Full_Parent_List.back()->parents()[0]);
        }
    }

        myfile << pid << ", " << prompt << ", " << !(apid == 12 || apid == 14 || apid == 16 || apid == 1000022 || apid == 1000039) << ", " <<  hp4.px() << ", " << hp4.py() << ", " << hp4.pz() << ", " << hp4.e() << ", " << abs(get_eta(gp)) << ",";
 
        for (auto event : Full_Parent_List)
        {
           myfile << event->pid() << ", " << event->momentum().e() << ", ";
        }
        myfile << "\n";*/
#endif



//////////////////////
    #ifdef ELECTRON_PARENTAGE_DEBUG
      if (apid == 11) print_history(gp);
    #endif

    // Store interacting prompt particles
    /// @todo Dress leptons?
    //Edited Tomek Procter Aug 2019 - Only stores prompt leptons.
    if ((apid == 22 || apid == 11 || apid == 13 || apid == 15) && prompt)
    {
      Particle* p = new Particle(p4, pid); // the event will take ownership of this pointer
      p->set_prompt(true);
      evt.add_particle(p);
      #ifdef JET_CHECKING
      if (apid == 22)
      {
        std::cout << "Prompt Photon: pz = " << p4.pz() << std::endl;
      }
      #endif
    }
    // Aggregate missing ET
    else if (apid == 12 || apid == 14 || apid == 16 || apid == 1000022 || apid == 1000039)
    {
      Particle* p = new Particle(p4, pid); // the event will take ownership of this pointer
      p->set_prompt(true);
      evt.add_particle(p);
      #ifdef MET_DEBUG
      std::cout << "Adding particle to met - PID: " << pid << "; Pz: " << p4.pz() << ", m: " << p4.m() <<std::endl;
      #endif
      vmet += p4;
    }
    // Store non-prompt momenta for Jet building
    if ((apid != MCUtils::PID::MUON) && !(apid == 12 || apid == 14 || apid == 16 || apid == 1000022 || apid == 1000039)) //Tomek Procter August 2019 - Muons don't go into jets, right? So what happens to non-prompt muons?
    {//Now added exclusion of invisibles as n longer else if as we want prompt gammas/leptons.
      PseudoJet pj = PseudoJet(hp4.px(), hp4.py(), hp4.pz(), hp4.e());
      //pj.set_user_index(apid);
      jetparticles.push_back(pj);
#ifdef JET_CHECKING
    if (TomeksCounter2 == 22)
    {
      jetfile << pid << ", " << pj.px() << ", " << pj.py() << ", " << pj.pz() << ", " << pj.e() << ", " << hp4.e() << ", " << pj.rap() << ",\n";
    }
#endif      
    }
    else
    {
       //cout << "LEFT OVER: " << apid << endl;
    }
  }
#ifdef SINGLE_EVENT_DEBUG
  myfile.close();
#endif
#ifdef JET_CHECKING
if (TomeksCounter2 == 22)
  {jetfile.close();}
#endif

  // MET
  #ifdef MET_DEBUG
  std::cout << "FINAL MISSING MOMENTUM IS: (" << vmet.px() <<", "<<vmet.py()<<", "<<vmet.pz()<<", "<<vmet.m()<<")"<<std::endl;
  #endif
  vmet.setPz(0);
  vmet.setM(0);
  #ifdef MET_DEBUG
  std::cout << "FINAL MISSING MOMENTUM IS: (" << vmet.px() <<", "<<vmet.py()<<", "<<vmet.pz()<<", "<<vmet.m()<<")"<<std::endl;
  #endif
  evt.set_missingmom(vmet);
  

  // Jets
  //Sorting the vector of pseudojets by z momentum so gambit and CBS have the same event order - makes debugging easier.
  //std::sort(jetparticles.begin(), jetparticles.end(), compare_particles_by_pz);
  //int TP_TEMP_COUNTER = 0;
  //for (auto particle : jetparticles)
  //{
    //std::cout << "Particle Number: " << TP_TEMP_COUNTER++ << "; Px: " << particle.px() << "; Pz: " << particle.pz() << "; e: " << particle.e() << std::endl;
  //}
  

  vector<PseudoJet> jets = get_jets(jetparticles, 0.4, 10.0); //Cut off momentum edited to match gambit at 10.0.
  

  for (const PseudoJet& pj : jets)
  {
    bool hasC = false, hasB = false, hastau = false;
    HEPUtils::P4 jetMom = HEPUtils::mk_p4(pj);
    //Tomek Procter Aug 2019 - replacing the block below with the logic from Gambit.
    /// @todo Bug in HEPUtils::get_jets means that constituent info is lost for now...
    // for (const PseudoJet& c : pj.constituents()) {
    //   if (c.user_index() == 4) hasC = true;
    //   if (c.user_index() == 5) hasB = true;
    // }
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

#ifdef NJET_DATA_INPUT
if (evt.jets().size() != gambit_njets[TomeksCounter2])
{
  std::cout << "Event " << TomeksCounter2 << " has jet discrepancy!" <<std::endl;
  std::cout << "Gambit Events: " << gambit_njets[TomeksCounter2]  << ", CBS njets: " << evt.jets().size() <<std::endl;
}
TomeksCounter2++;
#endif

#ifdef NJET_DATA_OUTPUT
    double sumpT = 0;
    for (int i{0}; i < evt.jets().size(); i++)
    {
      sumpT += evt.jets()[i]->pT();
    }
    std::fstream njetdatafile;
    njetdatafile.open("CBSNjetFile.csv", std::fstream::app);
    njetdatafile << evt.jets().size() << ", " << evt.met() << ", " << evt.jets()[0]->pT() << ", " << sumpT << "\n";
    njetdatafile.close();
#endif

#ifdef SINGLE_EVENT_DEBUG
  ++TomeksCounter2;
  //#ifdef COLLIDERBIT_DEBUG TP AUG 2019 - Lets get as much info out as possible!!!
    // Print event summary
    if (TomeksCounter2 >= 500)
    {
      cout << TomeksCounter2 << ". CBS Event Information\n";
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
    }
#endif

}

#endif
