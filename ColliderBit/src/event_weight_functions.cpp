//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  ColliderBit module functions for calculating 
///  event weights
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

    /// Set event weight to unity
    void _setEventWeightToUnity(HEPUtils::Event& event, double w)
    {
      std::cout << DEBUG_PREFIX << "This is _setEventWeightToUnity: setting event weight to " << w << std::endl;
      event.set_weight(w);
    }

    /// Module function providing an instance of EventWeightFunctionType
    /// pointing to _setEventWeightToUnity
    void setEventWeightToUnity(EventWeightFunctionType& result)
    {
      using namespace Pipes::setEventWeightToUnity;

      result = std::bind(_setEventWeightToUnity, std::placeholders::_1, 1.0);
    }


    /// Set event weight by cross-section
    void _setEventWeightByCrossSection(HEPUtils::Event& event, double w)
    {
      std::cout << DEBUG_PREFIX << "This is _setEventWeightByCrossSection: setting event weight to " << w << std::endl;
      event.set_weight(w);
    }

    /// Module function providing an instance of EventWeightFunctionType
    /// pointing to _setEventWeightByCrossSection
    void setEventWeightByCrossSection(EventWeightFunctionType& result)
    {
      using namespace Pipes::setEventWeightByCrossSection;
      std::cout << DEBUG_PREFIX << "This is setEventWeightByCrossSection: setting result to point to_setEventWeightByCrossSection " << std::endl;
      result = std::bind(_setEventWeightByCrossSection, std::placeholders::_1, 2.0);
    }

  } 
} 


