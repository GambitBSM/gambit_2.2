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

//#define SINGLE_EVENT_DEBUG
//#define ELECTRON_PARENTAGE_DEBUG
//#define JET_CHECKING

//#define ELECTRON_PARENTAGE_DEBUG_TWO

//Tomek Procter Aug 19. Note the MCUtils isParton function only checks
// for quarks/gluons (while in pythia, the function used in Gambit inlcudes
//diquarks too), so I'm manually defining this function using the isParton
//is diquark options in MCUtils
inline bool isParton(int pid)
{
   return (MCUtils::PID::isParton(pid) || MCUtils::PID::isDiquark(pid));
}


//TP Aug 19 - Adapted from Pythia function of same name.
inline bool fromHadron(HepMC3::ConstGenParticlePtr gp, bool print_progress = false, int counter = 0)
{
      if (MCUtils::PID::isHadron(abs(gp->pid()))) return true;
      if (isParton(abs(gp->pid()))) return false; // stop the walking at the end of the hadron level
      auto parent_vector = (gp->parents());
      if (parent_vector.size() == 0) return false;
#ifdef ELECTRON_PARENTAGE_DEBUG_TWO
      if (print_progress && counter < 10)
      {
        std::cout << "\nCounter - " << counter << ": ";
        for (auto parent:parent_vector)
        {
          std::cout << "(" << parent->pid() << ", " << parent->momentum().e() << "), ";
        } 
      }
#endif
//      if (parent_vector.size() > 1) std::cout << "\n\n\aTwo Events in Parent vector\n\n";
      for (const HepMC3::ConstGenParticlePtr parent : parent_vector)
      {
        if (counter >= 100) return false;
        if (fromHadron(parent, print_progress, ++counter)) return true;
      }
      return false;
}

//Tomek Procter August 2019. Iterative is much cleaner than recursive (if it works).
inline bool fromHadron_iterative(HepMC3::ConstGenParticlePtr gp)
{
    if (MCUtils::PID::isHadron(abs(gp->pid()))) return true;
    if (isParton(abs(gp->pid()))) return false;
    if (gp->parent_event() == NULL) return false;
    auto parent = gp->parents()[0];

    int loop_counter{0};
    while (!(parent->parent_event() == NULL) && (loop_counter++ < 250))
    {
        parent = parent->parents()[0];
        if (MCUtils::PID::isHadron(abs(parent->pid()))) return true;
        if (isParton(abs(parent->pid()))) return false;
        if (loop_counter == 248) std::cout << "\a\nLoop Counter Exceeded!!!";
    }
    return false;
}

//Tomek Procter Aug 2019
//This is a "fake" function - doesn't actually calc anything, just
//mimics the results for fromHadron seen in gambit.
inline bool fromHadron_test(HepMC3::ConstGenParticlePtr gp)
{
    if (abs(gp->pid()) == 1000022) return false;
    else return true;
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
  P4 vmet;
  vector<PseudoJet> jetparticles;
  evt.set_weight(ge.weight());

#ifdef SINGLE_EVENT_DEBUG
   std::fstream myfile;
   myfile.open("CBS_SINGLE_EVENT_FILE.csv", std::fstream::app);
   myfile << "Particle ID, Prompt?, Visible?, Px, Py, Pz, Pe, Abs(Eta), Parents, Parent Pe, etc.\n";
#endif

#ifdef JET_CHECKING
    std::fstream jetfile;
    jetfile.open("CBS_JetFile.csv", std::fstream::app);
#endif

  // Loop over all particles in the event
  for (HepMC3::ConstGenParticlePtr gp : ge.particles())
  {
    // Get status and PID code
    const int st = gp->status();
    const int pid = gp->pid();
    const int apid = fabs(pid);

    // Use physical particles only
    if (st != 1) continue;

    // Get 4-momentum
    const HepMC3::FourVector& hp4 = gp->momentum();

    //We need to define p4 some other way - as mkXYZE doesn't work, for now I'll do
    //const P4 p4 = P4::mkXYZM(hp4.px(), hp4.py(), hp4.pz(), hp4.e());

    const P4 p4 = mk_p4(PseudoJet(hp4.px(), hp4.py(), hp4.pz(), hp4.e()));
    

    //Lets test out adding an eta cut, like they have in gambit:
    if (abs(get_eta(gp)) > 5.0)
    {
    vmet += p4;
    continue; 
    }
       


/////////////////////////////////



// Promptness: for leptons and photons we're only interested if they don't come from hadron/tau decays: Added Tomek Procter Aug 19
        const bool prompt = !fromHadron(gp, (apid == 11)); //&& !fromTau(i, pevt);
        //cout << "\nLINE: " << __LINE__;
        //const bool visible = MCUtils::PID::isStrongInteracting(apid) || MCUtils::PID::isEMInteracting(apid); //Not needed by CBS (yet)
        //cout << "\nLINE: " << __LINE__;

 
#ifdef SINGLE_EVENT_DEBUG
//Find the full parent vector:

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
        myfile << "\n";
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
    }
    // Aggregate missing ET
    else if (apid == 12 || apid == 14 || apid == 16 || apid == 1000022 || apid == 1000039)
    {
      Particle* p = new Particle(p4, pid); // the event will take ownership of this pointer
      p->set_prompt(true);
      evt.add_particle(p);
      vmet += p4;

    }
    // Store non-prompt momenta for Jet building
    else if (apid != MCUtils::PID::MUON) //Tomek Procter August 2019 - Muons don't go into jets, right? So what happens to non-prompt muons?
    {
      PseudoJet pj = PseudoJet(hp4.px(), hp4.py(), hp4.pz(), hp4.e());
      //pj.set_user_index(apid);
      jetparticles.push_back(pj);
#ifdef JET_CHECKING
      jetfile << pj.px() << ", " << pj.py() << ", " << pj.pz() << ", " << pj.e() << ", " << hp4.e() << ", " << pj.rap() << ",\n";
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
  jetfile.close();
#endif

  // MET
  vmet.setPz(0);
  vmet.setM(0);
  evt.set_missingmom(vmet);
  

  // Jets
  //Sorting the vector of pseudojets by z momentum so gambit and CBS have the same event order - makes debugging easier.
  //std::sort(jetparticles.begin(), jetparticles.end(), compare_particles_by_pz);
  /*int TP_TEMP_COUNTER = 0;
  for (auto particle : jetparticles)
    {
      std::cout << "Particle Number: " << TP_TEMP_COUNTER++ << "; Pz: " << particle.pz() << "; Px: " << particle.px() <<std::endl;
    }*/
  

  vector<PseudoJet> jets = get_jets(jetparticles, 0.4, 10.0); //Cut off momentum edited to match gambit at 10.0.
  

  for (const PseudoJet& pj : jets)
  {
    bool hasC = false, hasB = false;
    /// @todo Bug in HEPUtils::get_jets means that constituent info is lost for now...
    // for (const PseudoJet& c : pj.constituents()) {
    //   if (c.user_index() == 4) hasC = true;
    //   if (c.user_index() == 5) hasB = true;
    // }
    evt.add_jet(new Jet(mk_p4(pj), hasB, hasC));
  }


#ifdef SINGLE_EVENT_DEBUG
  //#ifdef COLLIDERBIT_DEBUG TP AUG 2019 - Lets get as much info out as possible!!!
    // Print event summary
    cout << " CBS Event Information\n";
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
  //#endif
#endif

}

#endif
