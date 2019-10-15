#ifndef __loaded_types_FlexibleSUSY_CMSSM_2_4_0_hpp__
#define __loaded_types_FlexibleSUSY_CMSSM_2_4_0_hpp__ 1

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#include "wrapper_Command_line_options.h"
#include "wrapper_QedQcd.h"
#include "wrapper_Spectrum_generator_problems.h"
#include "wrapper_Spectrum_generator_settings.h"
#include "identification.hpp"
#pragma GCC diagnostic pop

// Indicate which types are provided by this backend, and what the symbols of their factories are.
#define FlexibleSUSY_CMSSM_2_4_0_all_data \
  (( /*class*/(flexiblesusy)(Command_line_options),    /*constructors*/(("Factory_Command_line_options_0__BOSS_1",())) (("Factory_Command_line_options_1__BOSS_2",(int, char**))) )) \
  (( /*class*/(softsusy)(QedQcd),    /*constructors*/(("Factory_QedQcd_0__BOSS_3",())) )) \
  (( /*class*/(flexiblesusy)(Spectrum_generator_problems),    /*constructors*/(("Factory_Spectrum_generator_problems_0__BOSS_4",())) )) \
  (( /*class*/(flexiblesusy)(Spectrum_generator_settings),    /*constructors*/(("Factory_Spectrum_generator_settings_0__BOSS_5",())) )) \

// If the default version has been loaded, set it as default.
#if ALREADY_LOADED(CAT_3(BACKENDNAME,_,CAT(Default_,BACKENDNAME)))
  SET_DEFAULT_VERSION_FOR_LOADING_TYPES(BACKENDNAME,SAFE_VERSION,CAT(Default_,BACKENDNAME))
#endif

// Undefine macros to avoid conflict with other backends.
#include "gambit/Backends/backend_undefs.hpp"

#endif /* __loaded_types_FlexibleSUSY_CMSSM_2_4_0_hpp__ */
