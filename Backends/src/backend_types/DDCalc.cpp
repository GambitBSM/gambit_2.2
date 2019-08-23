//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  Helper types for DDCalc backend.
///
///  *********************************************
///
///  Authors (add name and date if you modify):
///
///  \author Ben Farmer
///          (b.farmer@imperial.ac.uk)
///  \date 2019 Jul
///
///  *************************

#include "gambit/Backends/backend_types/DDCalc.hpp"
#include "gambit/Utils/model_parameters.hpp"
#include "gambit/Utils/util_types.hpp" 
#include "gambit/Models/safe_param_map.hpp"
#include <sstream>

namespace Gambit
{
    /// Default NREO_DM_nucleon_couplings constructor
    NREO_DM_nucleon_couplings::NREO_DM_nucleon_couplings()
    {
       for(int i=1; i<=15; i++)
       {
          c0[i] = 0;
          c1[i] = 0;
       }
    }

    /// NREO_DM_nucleon_couplings constuctor from ModelParameters object
    NREO_DM_nucleon_couplings::NREO_DM_nucleon_couplings(const ModelParameters& pars)
    {
       for(int i=1; i<=15; i++)
       {
          std::stringstream ss0;
          std::stringstream ss1;
          ss0<<"c0_"<<i;
          ss1<<"c1_"<<i;
          if(not pars.has(ss0.str()) or not pars.has(ss1.str()))
          {
             std::stringstream msg;
             msg<<"Error constructing NREO_DM_nucleon_couplings from ModelParameters! The supplied ModelParameters object (model name: "<<pars.getModelName()<<") does not contain NREO coupling parameters! Please check that an NREO-related ModelParameters object has been supplied.";
             backend_error().raise(LOCAL_INFO, msg.str());
          }
          c0[i] = pars[ss0.str()];
          c1[i] = pars[ss1.str()];
       } 
    }

    /// NREO_DM_nucleon_couplings constuctor from functor 'Params', i.e. 'safe_param_map' used to hold collected model parameters 
    NREO_DM_nucleon_couplings::NREO_DM_nucleon_couplings(const Models::safe_param_map<safe_ptr<const double>>& pars)
    {
       for(int i=1; i<=15; i++)
       {
          std::stringstream ss0;
          std::stringstream ss1;
          ss0<<"c0_"<<i;
          ss1<<"c1_"<<i;
          if(pars.find(ss0.str())==pars.end() or pars.find(ss1.str())==pars.end())
          {
             std::stringstream msg;
             msg<<"Error constructing NREO_DM_nucleon_couplings from functor Params map! The supplied Params map does not contain NREO coupling parameters! Please check that an NREO-related model has been activated with ALLOW_MODELS in the rollcall declaration for this module function.";
             backend_error().raise(LOCAL_INFO, msg.str());
          }
          c0[i] = *pars[ss0.str()];
          c1[i] = *pars[ss1.str()];
       }  
    }

    /// Function to prettify retrieval of couplings (also helpful for looping over 1,0 isospin integers) 
    double NREO_DM_nucleon_couplings::c(int iso, int o) const
    {
       if(iso!=0 and iso!=1)
       {
          std::stringstream msg;
          msg<<"Invalid isospin index (first argument) received ("<<iso<<")! Isospin index must be either 0 or 1";
          backend_error().raise(LOCAL_INFO, msg.str());
       }

       if(o<1 or o>15)
       {
          std::stringstream msg;
          msg<<"Invalid NREO index (second argument) received ("<<o<<")! Operator index must be an integer in the range [1,15]";
          backend_error().raise(LOCAL_INFO, msg.str()); 
       }

       double result;
       if(iso==0)
       {
           result = c0.at(o);
       }
       else if(iso==1)
       {
           result = c1.at(o);
       }

       return result;
    }

}
