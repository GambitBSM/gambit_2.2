//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  A big map connecting final-state PID pairs 
///  settings for the Prospino backend
///  
///  *********************************************
///
///  Authors (add name if you modify):
///
///  \author Anders Kvellestad
///          (a.kvellestad@imperial.ac.uk)
///  \date 2019 Nov
///
///  *********************************************


#ifndef __PID_pairs_to_prospino_settings_hpp__
#define __PID_pairs_to_prospino_settings_hpp__

#include "gambit/ColliderBit/ColliderBit_types.hpp"

#pragma once

namespace Gambit
{
  namespace ColliderBit
  {

    // One giant map initializer:
    static const std::map<PID_pair, prospino_settings>  PID_pairs_to_prospino_settings {
      //
      // Prospino settings: inlo, isq_ng_in, icoll_in, energy_in, i_error_in, finalState, ipart1, ipart2, isquark1_in, isquark2_in
      // gg
      std::make_pair( PID_pair(1000021, 1000021), prospino_settings(1, 1, 1, 13000., 0, "gg", 1, 1, 0, 0) ),
      // sg
      std::make_pair( PID_pair(1000021, 1000004), prospino_settings(1, 1, 1, 13000., 0, "sg", 1, 1, 4, 4) ),
      std::make_pair( PID_pair(1000021, 1000003), prospino_settings(1, 1, 1, 13000., 0, "sg", 1, 1, 3, 3) ),
      std::make_pair( PID_pair(1000021, 1000001), prospino_settings(1, 1, 1, 13000., 0, "sg", 1, 1, 2, 2) ),
      std::make_pair( PID_pair(1000021, 1000002), prospino_settings(1, 1, 1, 13000., 0, "sg", 1, 1, 1, 1) ),
      std::make_pair( PID_pair(1000021, 2000002), prospino_settings(1, 1, 1, 13000., 0, "sg", 1, 1, 1, 1) ),
      std::make_pair( PID_pair(1000021, 2000001), prospino_settings(1, 1, 1, 13000., 0, "sg", 1, 1, 2, 2) ),
      std::make_pair( PID_pair(1000021, 2000003), prospino_settings(1, 1, 1, 13000., 0, "sg", 1, 1, 3, 3) ),
      std::make_pair( PID_pair(1000021, 2000004), prospino_settings(1, 1, 1, 13000., 0, "sg", 1, 1, 4, 4) ),
      // sb
      std::make_pair( PID_pair(1000004,-1000004), prospino_settings(1, 1, 1, 13000., 0, "sb", 1, 1, 4, 4) ),
      std::make_pair( PID_pair(1000004,-1000003), prospino_settings(1, 1, 1, 13000., 0, "sb", 1, 1, 4, 3) ),
      std::make_pair( PID_pair(1000004,-1000001), prospino_settings(1, 1, 1, 13000., 0, "sb", 1, 1, 4, 2) ),
      std::make_pair( PID_pair(1000004,-1000002), prospino_settings(1, 1, 1, 13000., 0, "sb", 1, 1, 4, 1) ),
      std::make_pair( PID_pair(1000004,-2000002), prospino_settings(1, 1, 1, 13000., 0, "sb", 1, 1, 4, 1) ),
      std::make_pair( PID_pair(1000004,-2000001), prospino_settings(1, 1, 1, 13000., 0, "sb", 1, 1, 4, 2) ),
      std::make_pair( PID_pair(1000004,-2000003), prospino_settings(1, 1, 1, 13000., 0, "sb", 1, 1, 4, 3) ),
      std::make_pair( PID_pair(1000004,-2000004), prospino_settings(1, 1, 1, 13000., 0, "sb", 1, 1, 4, 4) ),
      std::make_pair( PID_pair(1000003,-1000003), prospino_settings(1, 1, 1, 13000., 0, "sb", 1, 1, 3, 3) ),
      std::make_pair( PID_pair(1000003,-1000001), prospino_settings(1, 1, 1, 13000., 0, "sb", 1, 1, 3, 2) ),
      std::make_pair( PID_pair(1000003,-1000002), prospino_settings(1, 1, 1, 13000., 0, "sb", 1, 1, 3, 1) ),
      std::make_pair( PID_pair(1000003,-2000002), prospino_settings(1, 1, 1, 13000., 0, "sb", 1, 1, 3, 1) ),
      std::make_pair( PID_pair(1000003,-2000001), prospino_settings(1, 1, 1, 13000., 0, "sb", 1, 1, 3, 2) ),
      std::make_pair( PID_pair(1000003,-2000003), prospino_settings(1, 1, 1, 13000., 0, "sb", 1, 1, 3, 3) ),
      std::make_pair( PID_pair(1000003,-2000004), prospino_settings(1, 1, 1, 13000., 0, "sb", 1, 1, 3, 4) ),
      std::make_pair( PID_pair(1000001,-1000001), prospino_settings(1, 1, 1, 13000., 0, "sb", 1, 1, 2, 2) ),
      std::make_pair( PID_pair(1000001,-1000002), prospino_settings(1, 1, 1, 13000., 0, "sb", 1, 1, 2, 1) ),
      std::make_pair( PID_pair(1000001,-2000002), prospino_settings(1, 1, 1, 13000., 0, "sb", 1, 1, 2, 1) ),
      std::make_pair( PID_pair(1000001,-2000001), prospino_settings(1, 1, 1, 13000., 0, "sb", 1, 1, 2, 2) ),
      std::make_pair( PID_pair(1000001,-2000003), prospino_settings(1, 1, 1, 13000., 0, "sb", 1, 1, 2, 3) ),
      std::make_pair( PID_pair(1000001,-2000004), prospino_settings(1, 1, 1, 13000., 0, "sb", 1, 1, 2, 4) ),
      std::make_pair( PID_pair(1000002,-1000002), prospino_settings(1, 1, 1, 13000., 0, "sb", 1, 1, 1, 1) ),
      std::make_pair( PID_pair(1000002,-2000002), prospino_settings(1, 1, 1, 13000., 0, "sb", 1, 1, 1, 1) ),
      std::make_pair( PID_pair(1000002,-2000001), prospino_settings(1, 1, 1, 13000., 0, "sb", 1, 1, 1, 2) ),
      std::make_pair( PID_pair(1000002,-2000003), prospino_settings(1, 1, 1, 13000., 0, "sb", 1, 1, 1, 3) ),
      std::make_pair( PID_pair(1000002,-2000004), prospino_settings(1, 1, 1, 13000., 0, "sb", 1, 1, 1, 4) ),
      std::make_pair( PID_pair(2000002,-2000002), prospino_settings(1, 1, 1, 13000., 0, "sb", 1, 1, 1, 1) ),
      std::make_pair( PID_pair(2000002,-2000001), prospino_settings(1, 1, 1, 13000., 0, "sb", 1, 1, 1, 2) ),
      std::make_pair( PID_pair(2000002,-2000003), prospino_settings(1, 1, 1, 13000., 0, "sb", 1, 1, 1, 3) ),
      std::make_pair( PID_pair(2000002,-2000004), prospino_settings(1, 1, 1, 13000., 0, "sb", 1, 1, 1, 4) ),
      std::make_pair( PID_pair(2000001,-2000001), prospino_settings(1, 1, 1, 13000., 0, "sb", 1, 1, 2, 2) ),
      std::make_pair( PID_pair(2000001,-2000003), prospino_settings(1, 1, 1, 13000., 0, "sb", 1, 1, 2, 3) ),
      std::make_pair( PID_pair(2000001,-2000004), prospino_settings(1, 1, 1, 13000., 0, "sb", 1, 1, 2, 4) ),
      std::make_pair( PID_pair(2000003,-2000003), prospino_settings(1, 1, 1, 13000., 0, "sb", 1, 1, 3, 3) ),
      std::make_pair( PID_pair(2000003,-2000004), prospino_settings(1, 1, 1, 13000., 0, "sb", 1, 1, 3, 4) ),
      std::make_pair( PID_pair(2000004,-2000004), prospino_settings(1, 1, 1, 13000., 0, "sb", 1, 1, 4, 4) ),
      // ss
      std::make_pair( PID_pair(1000004, 1000004), prospino_settings(1, 1, 1, 13000., 0, "ss", 1, 1, 4, 4) ),
      std::make_pair( PID_pair(1000004, 1000003), prospino_settings(1, 1, 1, 13000., 0, "ss", 1, 1, 4, 3) ),
      std::make_pair( PID_pair(1000004, 1000001), prospino_settings(1, 1, 1, 13000., 0, "ss", 1, 1, 4, 2) ),
      std::make_pair( PID_pair(1000004, 1000002), prospino_settings(1, 1, 1, 13000., 0, "ss", 1, 1, 4, 1) ),
      std::make_pair( PID_pair(1000004, 2000002), prospino_settings(1, 1, 1, 13000., 0, "ss", 1, 1, 4, 1) ),
      std::make_pair( PID_pair(1000004, 2000001), prospino_settings(1, 1, 1, 13000., 0, "ss", 1, 1, 4, 2) ),
      std::make_pair( PID_pair(1000004, 2000003), prospino_settings(1, 1, 1, 13000., 0, "ss", 1, 1, 4, 3) ),
      std::make_pair( PID_pair(1000004, 2000004), prospino_settings(1, 1, 1, 13000., 0, "ss", 1, 1, 4, 4) ),
      std::make_pair( PID_pair(1000003, 1000003), prospino_settings(1, 1, 1, 13000., 0, "ss", 1, 1, 3, 3) ),
      std::make_pair( PID_pair(1000003, 1000001), prospino_settings(1, 1, 1, 13000., 0, "ss", 1, 1, 3, 2) ),
      std::make_pair( PID_pair(1000003, 1000002), prospino_settings(1, 1, 1, 13000., 0, "ss", 1, 1, 3, 1) ),
      std::make_pair( PID_pair(1000003, 2000002), prospino_settings(1, 1, 1, 13000., 0, "ss", 1, 1, 3, 1) ),
      std::make_pair( PID_pair(1000003, 2000001), prospino_settings(1, 1, 1, 13000., 0, "ss", 1, 1, 3, 2) ),
      std::make_pair( PID_pair(1000003, 2000003), prospino_settings(1, 1, 1, 13000., 0, "ss", 1, 1, 3, 3) ),
      std::make_pair( PID_pair(1000003, 2000004), prospino_settings(1, 1, 1, 13000., 0, "ss", 1, 1, 3, 4) ),
      std::make_pair( PID_pair(1000001, 1000001), prospino_settings(1, 1, 1, 13000., 0, "ss", 1, 1, 2, 2) ),
      std::make_pair( PID_pair(1000001, 1000002), prospino_settings(1, 1, 1, 13000., 0, "ss", 1, 1, 2, 1) ),
      std::make_pair( PID_pair(1000001, 2000002), prospino_settings(1, 1, 1, 13000., 0, "ss", 1, 1, 2, 1) ),
      std::make_pair( PID_pair(1000001, 2000001), prospino_settings(1, 1, 1, 13000., 0, "ss", 1, 1, 2, 2) ),
      std::make_pair( PID_pair(1000001, 2000003), prospino_settings(1, 1, 1, 13000., 0, "ss", 1, 1, 2, 3) ),
      std::make_pair( PID_pair(1000001, 2000004), prospino_settings(1, 1, 1, 13000., 0, "ss", 1, 1, 2, 4) ),
      std::make_pair( PID_pair(1000002, 1000002), prospino_settings(1, 1, 1, 13000., 0, "ss", 1, 1, 1, 1) ),
      std::make_pair( PID_pair(1000002, 2000002), prospino_settings(1, 1, 1, 13000., 0, "ss", 1, 1, 1, 1) ),
      std::make_pair( PID_pair(1000002, 2000001), prospino_settings(1, 1, 1, 13000., 0, "ss", 1, 1, 1, 2) ),
      std::make_pair( PID_pair(1000002, 2000003), prospino_settings(1, 1, 1, 13000., 0, "ss", 1, 1, 1, 3) ),
      std::make_pair( PID_pair(1000002, 2000004), prospino_settings(1, 1, 1, 13000., 0, "ss", 1, 1, 1, 4) ),
      std::make_pair( PID_pair(2000002, 2000002), prospino_settings(1, 1, 1, 13000., 0, "ss", 1, 1, 1, 1) ),
      std::make_pair( PID_pair(2000002, 2000001), prospino_settings(1, 1, 1, 13000., 0, "ss", 1, 1, 1, 2) ),
      std::make_pair( PID_pair(2000002, 2000003), prospino_settings(1, 1, 1, 13000., 0, "ss", 1, 1, 1, 3) ),
      std::make_pair( PID_pair(2000002, 2000004), prospino_settings(1, 1, 1, 13000., 0, "ss", 1, 1, 1, 4) ),
      std::make_pair( PID_pair(2000001, 2000001), prospino_settings(1, 1, 1, 13000., 0, "ss", 1, 1, 2, 2) ),
      std::make_pair( PID_pair(2000001, 2000003), prospino_settings(1, 1, 1, 13000., 0, "ss", 1, 1, 2, 3) ),
      std::make_pair( PID_pair(2000001, 2000004), prospino_settings(1, 1, 1, 13000., 0, "ss", 1, 1, 2, 4) ),
      std::make_pair( PID_pair(2000003, 2000003), prospino_settings(1, 1, 1, 13000., 0, "ss", 1, 1, 3, 3) ),
      std::make_pair( PID_pair(2000003, 2000004), prospino_settings(1, 1, 1, 13000., 0, "ss", 1, 1, 3, 4) ),
      std::make_pair( PID_pair(2000004, 2000004), prospino_settings(1, 1, 1, 13000., 0, "ss", 1, 1, 4, 4) ),
      // tb
      std::make_pair( PID_pair(1000006,-1000006), prospino_settings(1, 1, 1, 13000., 0, "tb", 1, 1, 0, 0) ),
      std::make_pair( PID_pair(2000006,-2000006), prospino_settings(1, 1, 1, 13000., 0, "tb", 2, 2, 0, 0) ),
      // bb
      std::make_pair( PID_pair(1000005,-1000005), prospino_settings(1, 1, 1, 13000., 0, "bb", 1, 1, 0, 0) ),
      std::make_pair( PID_pair(2000005,-2000005), prospino_settings(1, 1, 1, 13000., 0, "bb", 2, 2, 0, 0) ),
      // nn
      std::make_pair( PID_pair(1000022, 1000022), prospino_settings(1, 1, 1, 13000., 0, "nn", 1, 1, 0, 0) ),
      std::make_pair( PID_pair(1000022, 1000023), prospino_settings(1, 1, 1, 13000., 0, "nn", 1, 2, 0, 0) ),
      std::make_pair( PID_pair(1000022, 1000025), prospino_settings(1, 1, 1, 13000., 0, "nn", 1, 3, 0, 0) ),
      std::make_pair( PID_pair(1000022, 1000035), prospino_settings(1, 1, 1, 13000., 0, "nn", 1, 4, 0, 0) ),
      std::make_pair( PID_pair(1000022, 1000024), prospino_settings(1, 1, 1, 13000., 0, "nn", 1, 5, 0, 0) ),
      std::make_pair( PID_pair(1000022, 1000037), prospino_settings(1, 1, 1, 13000., 0, "nn", 1, 6, 0, 0) ),
      std::make_pair( PID_pair(1000022,-1000024), prospino_settings(1, 1, 1, 13000., 0, "nn", 1, 7, 0, 0) ),
      std::make_pair( PID_pair(1000022,-1000037), prospino_settings(1, 1, 1, 13000., 0, "nn", 1, 8, 0, 0) ),
      std::make_pair( PID_pair(1000023, 1000023), prospino_settings(1, 1, 1, 13000., 0, "nn", 2, 2, 0, 0) ),
      std::make_pair( PID_pair(1000023, 1000025), prospino_settings(1, 1, 1, 13000., 0, "nn", 2, 3, 0, 0) ),
      std::make_pair( PID_pair(1000023, 1000035), prospino_settings(1, 1, 1, 13000., 0, "nn", 2, 4, 0, 0) ),
      std::make_pair( PID_pair(1000023, 1000024), prospino_settings(1, 1, 1, 13000., 0, "nn", 2, 5, 0, 0) ),
      std::make_pair( PID_pair(1000023, 1000037), prospino_settings(1, 1, 1, 13000., 0, "nn", 2, 6, 0, 0) ),
      std::make_pair( PID_pair(1000023,-1000024), prospino_settings(1, 1, 1, 13000., 0, "nn", 2, 7, 0, 0) ),
      std::make_pair( PID_pair(1000023,-1000037), prospino_settings(1, 1, 1, 13000., 0, "nn", 2, 8, 0, 0) ),
      std::make_pair( PID_pair(1000025, 1000025), prospino_settings(1, 1, 1, 13000., 0, "nn", 3, 3, 0, 0) ),
      std::make_pair( PID_pair(1000025, 1000035), prospino_settings(1, 1, 1, 13000., 0, "nn", 3, 4, 0, 0) ),
      std::make_pair( PID_pair(1000025, 1000024), prospino_settings(1, 1, 1, 13000., 0, "nn", 3, 5, 0, 0) ),
      std::make_pair( PID_pair(1000025, 1000037), prospino_settings(1, 1, 1, 13000., 0, "nn", 3, 6, 0, 0) ),
      std::make_pair( PID_pair(1000025,-1000024), prospino_settings(1, 1, 1, 13000., 0, "nn", 3, 7, 0, 0) ),
      std::make_pair( PID_pair(1000025,-1000037), prospino_settings(1, 1, 1, 13000., 0, "nn", 3, 8, 0, 0) ),
      std::make_pair( PID_pair(1000035, 1000035), prospino_settings(1, 1, 1, 13000., 0, "nn", 4, 4, 0, 0) ),
      std::make_pair( PID_pair(1000035, 1000024), prospino_settings(1, 1, 1, 13000., 0, "nn", 4, 5, 0, 0) ),
      std::make_pair( PID_pair(1000035, 1000037), prospino_settings(1, 1, 1, 13000., 0, "nn", 4, 6, 0, 0) ),
      std::make_pair( PID_pair(1000035,-1000024), prospino_settings(1, 1, 1, 13000., 0, "nn", 4, 7, 0, 0) ),
      std::make_pair( PID_pair(1000035,-1000037), prospino_settings(1, 1, 1, 13000., 0, "nn", 4, 8, 0, 0) ),
      std::make_pair( PID_pair(1000024, 1000024), prospino_settings(1, 1, 1, 13000., 0, "nn", 5, 5, 0, 0) ),
      std::make_pair( PID_pair(1000024, 1000037), prospino_settings(1, 1, 1, 13000., 0, "nn", 5, 6, 0, 0) ),
      std::make_pair( PID_pair(1000024,-1000024), prospino_settings(1, 1, 1, 13000., 0, "nn", 5, 7, 0, 0) ),
      std::make_pair( PID_pair(1000024,-1000037), prospino_settings(1, 1, 1, 13000., 0, "nn", 5, 8, 0, 0) ),
      std::make_pair( PID_pair(1000037, 1000037), prospino_settings(1, 1, 1, 13000., 0, "nn", 6, 6, 0, 0) ),
      std::make_pair( PID_pair(1000037,-1000024), prospino_settings(1, 1, 1, 13000., 0, "nn", 6, 7, 0, 0) ),
      std::make_pair( PID_pair(1000037,-1000037), prospino_settings(1, 1, 1, 13000., 0, "nn", 6, 8, 0, 0) ),
      // ll
      std::make_pair( PID_pair(1000011,-1000011), prospino_settings(1, 1, 1, 13000., 0, "ll", 1, 1, 0, 0) ),
      std::make_pair( PID_pair(2000011,-2000011), prospino_settings(1, 1, 1, 13000., 0, "ll", 2, 1, 0, 0) ),
      std::make_pair( PID_pair(1000012,-1000012), prospino_settings(1, 1, 1, 13000., 0, "ll", 3, 1, 0, 0) ),
      std::make_pair( PID_pair(1000012,-1000011), prospino_settings(1, 1, 1, 13000., 0, "ll", 4, 1, 0, 0) ),
      std::make_pair( PID_pair(1000011,-1000012), prospino_settings(1, 1, 1, 13000., 0, "ll", 5, 1, 0, 0) ),
      std::make_pair( PID_pair(1000015,-1000015), prospino_settings(1, 1, 1, 13000., 0, "ll", 6, 1, 0, 0) ),
      std::make_pair( PID_pair(2000015,-2000015), prospino_settings(1, 1, 1, 13000., 0, "ll", 7, 1, 0, 0) ),
      std::make_pair( PID_pair(1000015,-2000015), prospino_settings(1, 1, 1, 13000., 0, "ll", 8, 1, 0, 0) ),
      std::make_pair( PID_pair(1000016,-1000016), prospino_settings(1, 1, 1, 13000., 0, "ll", 9, 1, 0, 0) ),
      std::make_pair( PID_pair(1000016,-1000015), prospino_settings(1, 1, 1, 13000., 0, "ll", 10, 1, 0, 0) ),
      std::make_pair( PID_pair(1000015,-1000016), prospino_settings(1, 1, 1, 13000., 0, "ll", 11, 1, 0, 0) ),
      std::make_pair( PID_pair(1000016,-2000015), prospino_settings(1, 1, 1, 13000., 0, "ll", 12, 1, 0, 0) ),
      std::make_pair( PID_pair(2000015,-1000016), prospino_settings(1, 1, 1, 13000., 0, "ll", 13, 1, 0, 0) ),
    };

  }
}



#endif /* defined __PID_pairs_to_prospino_settings_hpp__ */
