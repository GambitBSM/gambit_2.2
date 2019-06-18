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


namespace Gambit
{

  namespace ColliderBit
  {

    /// Retrieve a BuckFast sim of ATLAS
    void getBuckFastATLAS(BaseDetector* &result)
    {
      using namespace Pipes::getBuckFastATLAS;
      static std::unique_ptr<BuckFast[]> buckies(new BuckFast[omp_get_max_threads()]);
      if (*Loop::iteration == START_SUBPROCESS)
      {
        BuckFast* bucky = &buckies[omp_get_thread_num()];
        bucky->smearElectronEnergy = &ATLAS::smearElectronEnergy;
        bucky->smearMuonMomentum   = &ATLAS::smearMuonMomentum;
        bucky->smearTaus           = &ATLAS::smearTaus;
        bucky->smearJets           = &ATLAS::smearJets;
        result = bucky;
      }
    }

    /// Retrieve a BuckFast sim of CMS
    void getBuckFastCMS(BaseDetector* &result)
    {
      using namespace Pipes::getBuckFastCMS;
      static std::unique_ptr<BuckFast[]> buckies(new BuckFast[omp_get_max_threads()]);
      if (*Loop::iteration == START_SUBPROCESS)
      {
        BuckFast* bucky = &buckies[omp_get_thread_num()];
        bucky->smearElectronEnergy = &CMS::smearElectronEnergy;
        bucky->smearMuonMomentum   = &CMS::smearMuonMomentum;
        bucky->smearTaus           = &CMS::smearTaus;
        bucky->smearJets           = &CMS::smearJets;
        result = bucky;
      }
    }

    /// Retrieve an Identity BuckFast sim (no sim)
    void getBuckFastIdentity(BaseDetector* &result)
    {
      using namespace Pipes::getBuckFastIdentity;
      static std::unique_ptr<BuckFast[]> buckies(new BuckFast[omp_get_max_threads()]);
      if (*Loop::iteration == START_SUBPROCESS) result = &buckies[omp_get_thread_num()];
    }

  }

}
