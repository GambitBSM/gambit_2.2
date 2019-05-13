// -*- C++ -*-
//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///
///  lhef2heputils: a Les Houches Event Format (LHEF)
///  -> HEPUtils::Event MC generator event file
///  converter, based on lhef2hepmc.
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
///  Hat tip: Leif Lonnblad for writing the LHEF
///  parser that actually makes this possible!
///
///  *********************************************

#include "LHEF.h"
#include "HEPUtils/Event.h"

/// Extract an LHE event as a HEPUtils::Event
HEPUtils::Event get_HEPUtils_event(const LHEF::Reader&);