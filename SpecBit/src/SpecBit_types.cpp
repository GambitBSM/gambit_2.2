//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  Source code for types for module SpecBit.
///  For instructions on adding new types, see
///  the corresponding header.
///
///  *********************************************
///
///  Authors (add name and date if you modify):
///
///  \author Janina Renk
///          (janina.renk@fysik.su.se)
///  \date 2019 July
///  *********************************************

#include <string>
#include <iostream>
#include <stdlib.h>     /* malloc, free, rand */
#include <valarray>

#include "gambit/SpecBit/SpecBit_types.hpp"

namespace Gambit
{
  namespace SpecBit
  {

    // initialise all members to -1 => in case vevacious does not run
    // and/or the bounce actions are not calculated this will be printed
    // to the output file
    VevaciousResultContainer::VevaciousResultContainer():
      lifetime(-1),                thermalProbability(-1),
      bounceActionThreshold(-1),   bounceActionThresholdThermal(-1),
      bounceActionStraight(-1),    bounceActionStraightThermal(-1),
      firstPathFinder(-1),         firstPathFinderThermal(-1),
      secondPathFinder(-1),        secondPathFinderThermal(-1)
    { 

      std::cout << "Init Container, values are " << lifetime << " thermalProb "<< thermalProbability << std::endl;
    }

    void VevaciousResultContainer::reset_results()
    { 
      lifetime = -1;                
      thermalProbability = -1;
      bounceActionThreshold = -1;  
      bounceActionThresholdThermal = -1;
      bounceActionStraight = -1;   
      bounceActionStraightThermal = -1;
      firstPathFinder = -1;        
      firstPathFinderThermal = -1;
      secondPathFinder = -1;       
      secondPathFinderThermal = -1;
      vevaciousRunFlag = false;
    }

    VevaciousResultContainer::~VevaciousResultContainer()
    {
      std::cout << "AAHHHGHG I DIED" << std::endl;
    }

    void VevaciousResultContainer::set_bounceActionThreshold(double val, bool thermal)
    {
      if(thermal) {bounceActionThresholdThermal = val;}
      else        {bounceActionThreshold = val;}
    }
  
    void VevaciousResultContainer::set_bounceActionStraight(double val, bool thermal)
    {
      if(thermal) {bounceActionStraightThermal = val;}
      else        {bounceActionStraight = val;}
    }
  
    void VevaciousResultContainer::set_firstPathFinder(double val, bool thermal)
    {
      if(thermal) {firstPathFinderThermal = val;}
      else        {firstPathFinder = val;}
    }
    
    void VevaciousResultContainer::set_secondPathFinder(double val, bool thermal)
    {
      if(thermal) {secondPathFinderThermal = val;}
      else        {secondPathFinder = val;}
    }
    
    double VevaciousResultContainer::get_bounceActionThreshold(bool thermal)
    {
      if(thermal) {return bounceActionThresholdThermal;}
      else        {return bounceActionThreshold;}
    }
  
    double VevaciousResultContainer::get_bounceActionStraight(bool thermal)
    {
      if(thermal) {return bounceActionStraightThermal;}
      else        {return bounceActionStraight;}
    }
  
    double VevaciousResultContainer::get_firstPathFinder(bool thermal)
    {
      if(thermal) {return firstPathFinderThermal;}
      else        {return firstPathFinder;}
    }
    
    double VevaciousResultContainer::get_secondPathFinder(bool thermal)
    {
      if(thermal) {return secondPathFinderThermal;}
      else        {return secondPathFinder;}
    }

  }
}
