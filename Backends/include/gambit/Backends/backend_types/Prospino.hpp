//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  Options struct for the Prospino backend
///
///  *********************************************
///
///  Authors (add name and date if you modify):
///  
///  \author Anders Kvellestad
///  \date 2019 Jun
///
///  *********************************************

#include <map>

#include "gambit/Utils/util_types.hpp"

#ifndef __PROSPINO_types_hpp__
#define __PROSPINO_types_hpp__

namespace Gambit
{

    struct prospino_settings
    {
        Finteger inlo;
        Finteger isq_ng_in;
        Finteger icoll_in;
        Fdouble energy_in;
        Finteger i_error_in;

        Fstring<2> final_state_in;
        Finteger ipart1_in;
        Finteger ipart2_in;
        Finteger isquark1_in;
        Finteger isquark2_in;
    };

}

#endif // defined __PROSPINO_types_hpp__
