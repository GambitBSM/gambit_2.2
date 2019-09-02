// -*- C++ -*-
//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///
///  lhef2heputils: a Les Houches Event Format (LHEF)
///  -> HEPUtils::Event MC generator event file
///  converter, based on lhef2hepmc.
///
///  Hat tip: Leif Lönnblad for writing the LHEF
///  parser that actually makes this possible!
///
///  *********************************************
///
///  Authors (add name and date if you modify):
///
///  \author Andy Buckley
///  (andy.buckley@cern.ch)
///  \date May 2019
///
///  \author Pat Scott
///  (p.scott@imperial.ac.uk)
///  \date May 2019
///
///  *********************************************

using namespace std;

#include <iostream>
#include "LHEF.h"
#include "HEPUtils/FastJet.h"
#include "gambit/ColliderBit/lhef2heputils.hpp"

//#define COLLIDERBIT_DEBUG

using namespace HEPUtils;
using namespace FJNS;

/// Extract an LHE event as a HEPUtils::Event
void get_HEPUtils_event(const LHEF::Reader& lhe, Event& evt)
{

  P4 vmet;
  vector<PseudoJet> jetparticles;

  evt.set_weight(lhe.hepeup.weight());

  // Loop over all particles in the event
  for (int i = 0; i < lhe.hepeup.NUP; ++i)
  {
    // Get status and PID code
    const int st = lhe.hepeup.ISTUP[i];
    const int apid = fabs(lhe.hepeup.IDUP[i]);

    // Use LHE-stable particles only
    if (st != 1) continue;

    // Get 4-momentum
    const P4 p4 = P4::mkXYZM(lhe.hepeup.PUP[i][0], lhe.hepeup.PUP[i][1], lhe.hepeup.PUP[i][2], lhe.hepeup.PUP[i][4]);

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