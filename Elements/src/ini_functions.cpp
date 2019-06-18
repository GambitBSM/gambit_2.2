//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  Functions for triggering initialisation code.
///
///  *********************************************
///
///  Authors
///  =======
///
///  (add name and date if you modify)
///
///  \author Pat Scott
///          (p.scott@imperial.ac.uk)
///  \date 2015 Feb
///
///  \author Peter Athron
///          (peter.athron@coepp.org.au)
///  \date 2015
///
///  \author Christoph Weniger
///          (c.weniger@uva.nl)
///  \date 2016 Feb
///
///  \author Tomas Gonzalo
///          (t.e.gonzalo@fys.uio.no)
///  \date 2016 Sep
///
///  *********************************************

#include "gambit/Elements/ini_functions.hpp"
#include "gambit/Elements/ini_catch.hpp"
#include "gambit/Elements/functors.hpp"
#include "gambit/Elements/equivalency_singleton.hpp"
#include "gambit/Backends/backend_singleton.hpp"
#include "gambit/Models/claw_singleton.hpp"
#include "gambit/Logs/logging.hpp"

namespace Gambit
{

  /// Helper function for adding a type equivalency at initialisation
  int add_equivrelation(str s1, str s2)
  {
    try
    {
      Utils::typeEquivalencies().add(s1,s2);
    }
    catch (std::exception& e) { ini_catch(e); }
    return 0;
  }

  /// Helper function for passing default backend information at initialisation
  int pass_default_to_backendinfo(str be, str def)
  {
    try
    {
      Backends::backendInfo().default_safe_versions[be] = def;
    }
    catch (std::exception& e) { ini_catch(e); }
    return 0;
  }

  /// Runtime addition of model to GAMBIT model database
  int add_model(str model, str parent)
  {
    try
    {
      Models::ModelDB().declare_model(model, parent);
    }
    catch (std::exception& e) { ini_catch(e); }
    return 0;
  }

  /// Add a new parameter to a primary model functor
  int add_parameter(model_functor& primary_parameters, str param)
  {
    try
    {
      primary_parameters.addParameter(param);
    }
    catch (std::exception& e) { ini_catch(e); }
    return 0;
  }

  /// Set model name string in a primary model functor
  int set_model_name(model_functor& primary_parameters, str model_name)
  {
    try
    {
      primary_parameters.setModelName(model_name);
    }
    catch (std::exception& e) { ini_catch(e); }
    return 0;
  }

  /// Tell a model functor to take its parameter definition from another model functor.
  int copy_parameters(model_functor& donor, model_functor& donee, bool add_friend, str model, str model_x)
  {
    try
    {
      donor.donateParameters(donee);
      if (add_friend) Models::ModelDB().add_friend(model, model_x);
    }
    catch (std::exception& e) { ini_catch(e); }
    return 0;
  }

  /// Register a model functor.
  int register_model_functor(std::map<str, bool(*)()> map_bools, std::map<str, str> iCanDo,
   bool(*provides_function)(), str function, str capability)
  {
    try
    {
      map_bools[capability] = provides_function;
      iCanDo[function] = "ModelParameters";
    }
    catch (std::exception& e) { ini_catch(e); }
    return 0;
  }

  /// Create a log tag for a new module.
  int register_module_with_log(str module)
  {
    int mytag;
    try
    {
      mytag = Logging::getfreetag();
      Logging::tag2str()[mytag] = module;
      Logging::components().insert(mytag);
    }
    catch (std::exception& e) { ini_catch(e); }
    return mytag;
  }

  /// Register a function with a module.
  int register_function(module_functor_common& f, bool can_manage, safe_ptr<bool>* done,
   std::map<str,str>& iCanDo, std::map<str, bool(*)()>& map, bool(&provides)(), safe_ptr<Options>& opts)
  {
    try
    {
      if (can_manage)
      {
        f.setCanBeLoopManager(true);
        *done = f.loopIsDone();
      }
      map[f.capability()] = &provides;
      iCanDo[f.name()] = f.type();
      opts = f.getOptions();
    }
    catch (std::exception& e) { ini_catch(e); }
    return 0;
  }

}
