//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  ColliderBit event loop functions returning
///  detector simulations.
///
///  *********************************************
///
///  Authors (add name and date if you modify):
///
///  \author Abram Krislock
///          (a.m.b.krislock@fys.uio.no)
///
///  \author Aldo Saavedra
///
///  \author Andy Buckley
///
///  \author Chris Rogan
///          (crogan@cern.ch)
///  \date 2014 Aug
///  \date 2015 May
///
///  \author Pat Scott
///          (p.scott@imperial.ac.uk)
///  \date 2015 Jul
///  \date 2018 Jan
///  \date 2019 Jan
///
///  \author Anders Kvellestad
///          (anders.kvellestad@fys.uio.no)
///  \date   2017 March
///  \date   2018 Jan
///  \date   2018 May
///
///  *********************************************

#include <memory>

#include "gambit/ColliderBit/ColliderBit_eventloop.hpp"
#include "gambit/ColliderBit/ATLASEfficiencies.hpp"
#include "gambit/ColliderBit/CMSEfficiencies.hpp"
#include "gambit/ColliderBit/detectors/BuckFast.hpp"

#include "HEPUtils/FastJet.h"

// #define COLLIDERBIT_DEBUG

namespace Gambit
{

  namespace ColliderBit
  {

    /// Get a BuckFast detector simulation
    BaseDetector* getBuckFast(const str& detname
                              ,int iteration
                              //,const Options& runOptions
                              )
    {
      // Where's my Bucky?
      static std::unique_ptr<BuckFast[]> bucky(new BuckFast[omp_get_max_threads()]);
      int mine = omp_get_thread_num();

      if (iteration == START_SUBPROCESS)
      {
        // Assign detector functions
        if (detname == "ATLAS")
        {
          bucky[mine].smearElectronEnergy = &ATLAS::smearElectronEnergy;
          bucky[mine].smearMuonMomentum   = &ATLAS::smearMuonMomentum ;
          bucky[mine].smearTaus           = &ATLAS::smearTaus;
          bucky[mine].smearJets           = &ATLAS::smearJets;
        }
        else if (detname == "CMS")
        {
          bucky[mine].smearElectronEnergy = &CMS::smearElectronEnergy;
          bucky[mine].smearMuonMomentum   = &CMS::smearMuonMomentum;
          bucky[mine].smearTaus           = &CMS::smearTaus;
          bucky[mine].smearJets           = &CMS::smearJets;
        }
        else if (detname == "Identity") { /* relax */ }
        else
        {
          ColliderBit_error().raise(LOCAL_INFO, "Unrecognised detector name.");
        }
      }

      // Paper-bag it
      return &bucky[mine];

    }

    /// Retrieve a BuckFast sim of EXPERIMENT
    #define GET_BUCKFAST_AS_BASE_DETECTOR(NAME, EXPERIMENT)             \
    void NAME(BaseDetector* &result)                                    \
    {                                                                   \
      using namespace Pipes::NAME;                                      \
      result = getBuckFast(#EXPERIMENT, *Loop::iteration/*, *runOptions*/); \
    }

    GET_BUCKFAST_AS_BASE_DETECTOR(getBuckFastATLAS, ATLAS)
    GET_BUCKFAST_AS_BASE_DETECTOR(getBuckFastCMS, CMS)
    GET_BUCKFAST_AS_BASE_DETECTOR(getBuckFastIdentity, Identity)

  }

}
