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
///  *********************************************

#include "gambit/cmake/cmake_variables.hpp"

#ifndef EXCLUDE_HEPMC

using namespace std;

#include <iostream>

#include "gambit/ColliderBit/hepmc2heputils.hpp"

#include "HepMC3/GenEvent.h"
#include "HepMC3/GenParticle.h"
#include "HEPUtils/FastJet.h"

//#define COLLIDERBIT_DEBUG

using namespace HEPUtils;
using namespace FJNS;

/// Extract a HepMC event as a HEPUtils::Event
void get_HEPUtils_event(const HepMC3::GenEvent& ge, HEPUtils::Event& evt)
{

  P4 vmet;
  vector<PseudoJet> jetparticles;

  evt.set_weight(ge.weight());

  // Loop over all particles in the event
  for (HepMC3::ConstGenParticlePtr gp : ge.particles())
  {
    // Get status and PID code
    const int st = gp->status();
    const int apid = fabs(gp->pid());

    // Use physical particles only
    if (st != 1 && st != 2) continue;

    // Get 4-momentum
    const HepMC3::FourVector& hp4 = gp->momentum();
    const P4 p4 = P4::mkXYZM(hp4.px(), hp4.py(), hp4.pz(), hp4.e());

    // Store interacting prompt particles
    /// @todo Dress leptons?
    if (apid == 22 || apid == 11 || apid == 13 || apid == 15)
    {
      evt.add_particle(new Particle(p4, apid));
    }

    // Aggregate missing ET
    else if (apid == 12 || apid == 14 || apid == 16 || apid == 1000022 || apid == 1000039)
    {
      evt.add_particle(new Particle(p4, apid));
      vmet += p4;
    }

    // Store non-prompt momenta for Jet building
    else
    {
      PseudoJet pj = mk_pj(p4);
      pj.set_user_index(apid);
      jetparticles.push_back(pj);
    }
  }

  // MET
  vmet.setPz(0);
  vmet.setM(0);
  evt.set_missingmom(vmet);

  // Jets
  vector<PseudoJet> jets = get_jets(jetparticles, 0.4, 20.0);
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

  #ifdef COLLIDERBIT_DEBUG
    // Print event summary
    cout << "  MET  = " << evt.met() << " GeV" << endl;
    cout << "  #e   = " << evt.electrons().size() << endl;
    cout << "  #mu  = " << evt.muons().size() << endl;
    cout << "  #tau = " << evt.taus().size() << endl;
    cout << "  #jet = " << evt.jets().size() << endl;
    cout << "  #pho  = " << evt.photons().size() << endl;
    cout << endl;
  #endif

}

#endif
