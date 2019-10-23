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

#include "gambit/cmake/cmake_variables.hpp"

#ifndef EXCLUDE_HEPMC

#include "HEPUtils/Event.h"

/// Forward declaration to cut down on includes
namespace HepMC3
{
  class GenEvent;
}

/// Extract a HepMC event as a HEPUtils::Event
void get_HEPUtils_event(const HepMC3::GenEvent&, HEPUtils::Event&, double);

#endif
