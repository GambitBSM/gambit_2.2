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
///  \date May 2019
///
///  *********************************************

#include "HEPUtils/Event.h"

/// Forward declaration to cut down on includes
namespace HepMC {
  class GenEvent;
  class IO_GenEvent;
}

/// Extract a HepMC event as a HEPUtils::Event
void get_HEPUtils_event(std::unique_ptr<const HepMC::GenEvent>, HEPUtils::Event&);
