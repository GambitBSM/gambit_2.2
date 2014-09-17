//  *********************************************
//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  Macros for declaring different types for 
///  GAMBIT.  Version to be included in all
///  compilation units other than the main 
///  program.
///
///  *********************************************
///
///  Authors (add name and date if you modify):
///   
///  \author Pat Scott  
///          (p.scott@imperial.ac.uk)
///  \date 2014 Sep
///
///  *********************************************

#ifndef __type_macros_hpp__
#define __type_macros_hpp__

#include "util_macros.hpp"
#include <boost/preprocessor/seq/for_each.hpp>

/// Set default backend version for BOSSed types.
#define MAKE_DEFAULT_VERSION_FOR_LOADING_TYPES(BE,VER)                              \
 BOOST_PP_SEQ_FOR_EACH(TYPEDEFAULT, CAT_3(BE,_,VER), CAT_5(BE,_,VER,_,all_types))       

/// Helper macro for setting default backend version for BOSSed types.
#define TYPEDEFAULT(r,NSPACE,TNAME) namespace Gambit { using NSPACE::TNAME; }                                           


#endif //__type_macros_hpp__


