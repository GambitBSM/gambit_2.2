//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  ColliderBit module functions for calculating 
///  event weights
///
///  The weight functions in this file are
///  independent of the particular Py8Collider type
///  used in event generation.
///
///  *********************************************
///
///  Authors (add name and date if you modify):
///
///  \author Anders Kvellestad
///          (a.kvellestad@imperial.ac.uk)
///  \date   2019 Sept
///
///  *********************************************

#include "gambit/ColliderBit/ColliderBit_eventloop.hpp"

#define COLLIDERBIT_DEBUG
#define DEBUG_PREFIX "DEBUG: OMP thread " << omp_get_thread_num() << ":  "

namespace Gambit
{

  namespace ColliderBit
  {

    /// A function that sets the event weight to unity
    void _setEventWeight_unity(HEPUtils::Event& event, const BaseCollider*)  // <-- Ignoring second argument
    {
      event.set_weight(1.0);
    }

    /// Module function providing an instance of EventWeighterType_Py8Collider
    /// pointing to _setEventWeight_unity
    void setEventWeight_unity(EventWeighterType_Py8Collider& result)
    {
      using namespace Pipes::setEventWeight_unity;

      result = std::bind(_setEventWeight_unity, std::placeholders::_1, std::placeholders::_2);
    }


  } 
} 


