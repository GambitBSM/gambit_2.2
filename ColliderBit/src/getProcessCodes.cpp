//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  ColliderBit module functions for dealing with
///  lists of the active processes
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


    // enum specialIterations { BASE_INIT = -1,
    //                          COLLIDER_INIT = -2,
    //                          START_SUBPROCESS = -3,
    //                          COLLECT_CONVERGENCE_DATA = -4,
    //                          CHECK_CONVERGENCE = -5,
    //                          END_SUBPROCESS = -6,
    //                          COLLIDER_FINALIZE = -7,
    //                          BASE_FINALIZE = -8};


    /// Get the list of active Pythia process codes
    void getPythiaProcessCodes(std::vector<int>& result)
    {
      using namespace Pipes::getPythiaProcessCodes;

      if (*Loop::iteration == COLLIDER_INIT)
      {
        result.clear();
      }

      if (*Loop::iteration == XSEC_CALCULATION)
      {
        result = Dep::HardScatteringSim->codesHard();
      }

      // if (*Loop::iteration == BASE_INIT) { cout << DEBUG_PREFIX << "getPythiaProcessCodes: in BASE_INIT: result.size() = " << result.size() << endl; }
      // if (*Loop::iteration == COLLIDER_INIT) { cout << DEBUG_PREFIX << "getPythiaProcessCodes: in COLLIDER_INIT: result.size() = " << result.size() << endl; }
      // if (*Loop::iteration == COLLIDER_INIT_OMP) { cout << DEBUG_PREFIX << "getPythiaProcessCodes: in COLLIDER_INIT_OMP: result.size() = " << result.size() << endl; }
      // if (*Loop::iteration == XSEC_CALCULATION) { cout << DEBUG_PREFIX << "getPythiaProcessCodes: in XSEC_CALCULATION: result.size() = " << result.size() << endl; }
      // if (*Loop::iteration == START_SUBPROCESS) { cout << DEBUG_PREFIX << "getPythiaProcessCodes: in START_SUBPROCESS: result.size() = " << result.size() << endl; }
      // if (*Loop::iteration == 0) { cout << DEBUG_PREFIX << "getPythiaProcessCodes: in ITERATION 0: result.size() = " << result.size() << endl; }
      // if (*Loop::iteration == 1) { cout << DEBUG_PREFIX << "getPythiaProcessCodes: in ITERATION 1: result.size() = " << result.size() << endl; }
      // if (*Loop::iteration == END_SUBPROCESS) { cout << DEBUG_PREFIX << "getPythiaProcessCodes: in END_SUBPROCESS: result.size() = " << result.size() << endl; }
      // if (*Loop::iteration == COLLIDER_FINALIZE) { cout << DEBUG_PREFIX << "getPythiaProcessCodes: in COLLIDER_FINALIZE: result.size() = " << result.size() << endl; }
      // if (*Loop::iteration == BASE_FINALIZE) { cout << DEBUG_PREFIX << "getPythiaProcessCodes: in BASE_FINALIZE: result.size() = " << result.size() << endl; }

    }


    /// Translate list of Pythia process codes to list of (PID,PID) pairs
    /// for the two final state particles of the hard process
    void getProcessPIDPairs(vec_PID_pairs& result)
    {
      using namespace Pipes::getProcessPIDPairs;

      static bool first = true;
      static std::multimap<int,PID_pair> process_to_PIDs;

      // Initialize the map from Pythia process codes to PID pairs the first time
      // this module function is run
      if (first)
      {
        // Gluino pair prodction
        process_to_PIDs.insert( std::pair<int,PID_pair>(1201, PID_pair(1000021, 1000021)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1202, PID_pair(1000021, 1000021)) );

        // Gluino--squark production
        process_to_PIDs.insert( std::pair<int,PID_pair>(1203, PID_pair(1000001, 1000021)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1203, PID_pair(-1000001, 1000021)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1204, PID_pair(1000002, 1000021)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1204, PID_pair(-1000002, 1000021)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1205, PID_pair(1000003, 1000021)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1205, PID_pair(-1000003, 1000021)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1206, PID_pair(1000004, 1000021)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1206, PID_pair(-1000004, 1000021)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1207, PID_pair(1000005, 1000021)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1207, PID_pair(-1000005, 1000021)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1209, PID_pair(2000001, 1000021)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1209, PID_pair(-2000001, 1000021)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1210, PID_pair(2000002, 1000021)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1210, PID_pair(-2000002, 1000021)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1211, PID_pair(2000003, 1000021)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1211, PID_pair(-2000003, 1000021)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1212, PID_pair(2000004, 1000021)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1212, PID_pair(-2000004, 1000021)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1213, PID_pair(2000005, 1000021)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1213, PID_pair(-2000005, 1000021)) );
        
        // Squark--anti-squark production
        process_to_PIDs.insert( std::pair<int,PID_pair>(1215, PID_pair(1000001, -1000001)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1216, PID_pair(1000002, -1000002)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1217, PID_pair(1000003, -1000003)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1218, PID_pair(1000004, -1000004)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1219, PID_pair(1000005, -1000005)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1220, PID_pair(1000006, -1000006)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1221, PID_pair(2000001, -2000001)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1222, PID_pair(2000002, -2000002)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1223, PID_pair(2000003, -2000003)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1224, PID_pair(2000004, -2000004)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1225, PID_pair(2000005, -2000005)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1226, PID_pair(2000006, -2000006)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1231, PID_pair(1000001, -1000001)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1232, PID_pair(1000001, -1000003)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1232, PID_pair(1000003, -1000001)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1233, PID_pair(1000003, -1000001)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1233, PID_pair(1000001, -1000003)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1234, PID_pair(1000001, -1000005)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1234, PID_pair(1000005, -1000001)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1235, PID_pair(1000005, -1000001)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1235, PID_pair(1000001, -1000005)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1236, PID_pair(1000001, -2000001)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1236, PID_pair(2000001, -1000001)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1237, PID_pair(2000001, -1000001)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1237, PID_pair(1000001, -2000001)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1238, PID_pair(1000001, -2000003)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1238, PID_pair(2000003, -1000001)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1239, PID_pair(2000003, -1000001)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1239, PID_pair(1000001, -2000003)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1240, PID_pair(1000001, -2000005)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1240, PID_pair(2000005, -1000001)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1241, PID_pair(2000005, -1000001)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1241, PID_pair(1000001, -2000005)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1242, PID_pair(1000002, -1000002)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1243, PID_pair(1000002, -1000004)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1243, PID_pair(1000004, -1000002)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1244, PID_pair(1000004, -1000002)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1244, PID_pair(1000002, -1000004)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1247, PID_pair(1000002, -2000002)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1247, PID_pair(2000002, -1000002)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1248, PID_pair(2000002, -1000002)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1248, PID_pair(1000002, -2000002)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1249, PID_pair(1000002, -2000004)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1249, PID_pair(2000004, -1000002)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1250, PID_pair(2000004, -1000002)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1250, PID_pair(1000002, -2000004)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1253, PID_pair(1000002, -1000001)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1253, PID_pair(-1000002, 1000001)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1254, PID_pair(1000002, -1000003)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1254, PID_pair(-1000002, 1000003)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1255, PID_pair(1000002, -1000005)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1255, PID_pair(-1000002, 1000005)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1256, PID_pair(1000002, -2000001)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1256, PID_pair(-1000002, 2000001)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1257, PID_pair(1000002, -2000003)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1257, PID_pair(-1000002, 2000003)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1258, PID_pair(1000002, -2000005)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1258, PID_pair(-1000002, 2000005)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1259, PID_pair(1000003, -1000003)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1260, PID_pair(1000003, -1000005)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1260, PID_pair(1000005, -1000003)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1261, PID_pair(1000005, -1000003)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1261, PID_pair(1000003, -1000005)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1262, PID_pair(1000003, -2000001)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1262, PID_pair(2000001, -1000003)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1263, PID_pair(2000001, -1000003)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1263, PID_pair(1000003, -2000001)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1264, PID_pair(1000003, -2000003)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1264, PID_pair(2000003, -1000003)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1265, PID_pair(2000003, -1000003)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1265, PID_pair(1000003, -2000003)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1266, PID_pair(1000003, -2000005)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1266, PID_pair(2000005, -1000003)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1267, PID_pair(2000005, -1000003)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1267, PID_pair(1000003, -2000005)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1268, PID_pair(1000004, -1000004)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1271, PID_pair(1000004, -2000002)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1271, PID_pair(2000002, -1000004)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1272, PID_pair(2000002, -1000004)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1272, PID_pair(1000004, -2000002)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1273, PID_pair(1000004, -2000004)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1273, PID_pair(2000004, -1000004)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1274, PID_pair(2000004, -1000004)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1274, PID_pair(1000004, -2000004)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1277, PID_pair(1000004, -1000001)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1277, PID_pair(-1000004, 1000001)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1278, PID_pair(1000004, -1000003)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1278, PID_pair(-1000004, 1000003)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1279, PID_pair(1000004, -1000005)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1279, PID_pair(-1000004, 1000005)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1280, PID_pair(1000004, -2000001)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1280, PID_pair(-1000004, 2000001)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1281, PID_pair(1000004, -2000003)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1281, PID_pair(-1000004, 2000003)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1282, PID_pair(1000004, -2000005)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1282, PID_pair(-1000004, 2000005)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1283, PID_pair(1000005, -1000005)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1284, PID_pair(1000005, -2000001)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1284, PID_pair(2000001, -1000005)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1285, PID_pair(2000001, -1000005)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1285, PID_pair(1000005, -2000001)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1286, PID_pair(1000005, -2000003)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1286, PID_pair(2000003, -1000005)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1287, PID_pair(2000003, -1000005)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1287, PID_pair(1000005, -2000003)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1288, PID_pair(1000005, -2000005)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1288, PID_pair(2000005, -1000005)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1289, PID_pair(2000005, -1000005)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1289, PID_pair(1000005, -2000005)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1290, PID_pair(1000006, -1000006)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1295, PID_pair(1000006, -2000006)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1295, PID_pair(2000006, -1000006)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1296, PID_pair(2000006, -1000006)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1296, PID_pair(1000006, -2000006)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1297, PID_pair(1000006, -1000001)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1297, PID_pair(-1000006, 1000001)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1298, PID_pair(1000006, -1000003)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1298, PID_pair(-1000006, 1000003)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1299, PID_pair(1000006, -1000005)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1299, PID_pair(-1000006, 1000005)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1302, PID_pair(1000006, -2000005)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1302, PID_pair(-1000006, 2000005)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1303, PID_pair(2000001, -2000001)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1304, PID_pair(2000001, -2000003)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1304, PID_pair(2000003, -2000001)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1305, PID_pair(2000003, -2000001)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1305, PID_pair(2000001, -2000003)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1306, PID_pair(2000001, -2000005)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1306, PID_pair(2000005, -2000001)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1307, PID_pair(2000005, -2000001)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1307, PID_pair(2000001, -2000005)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1308, PID_pair(2000002, -2000002)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1309, PID_pair(2000002, -2000004)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1309, PID_pair(2000004, -2000002)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1310, PID_pair(2000004, -2000002)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1310, PID_pair(2000002, -2000004)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1313, PID_pair(2000002, -1000001)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1313, PID_pair(-2000002, 1000001)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1314, PID_pair(2000002, -1000003)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1314, PID_pair(-2000002, 1000003)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1315, PID_pair(2000002, -1000005)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1315, PID_pair(-2000002, 1000005)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1316, PID_pair(2000002, -2000001)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1316, PID_pair(-2000002, 2000001)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1317, PID_pair(2000002, -2000003)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1317, PID_pair(-2000002, 2000003)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1318, PID_pair(2000002, -2000005)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1318, PID_pair(-2000002, 2000005)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1319, PID_pair(2000003, -2000003)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1320, PID_pair(2000003, -2000005)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1320, PID_pair(2000005, -2000003)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1321, PID_pair(2000005, -2000003)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1321, PID_pair(2000003, -2000005)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1322, PID_pair(2000004, -2000004)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1325, PID_pair(2000004, -1000001)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1325, PID_pair(-2000004, 1000001)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1326, PID_pair(2000004, -1000003)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1326, PID_pair(-2000004, 1000003)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1327, PID_pair(2000004, -1000005)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1327, PID_pair(-2000004, 1000005)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1328, PID_pair(2000004, -2000001)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1328, PID_pair(-2000004, 2000001)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1329, PID_pair(2000004, -2000003)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1329, PID_pair(-2000004, 2000003)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1330, PID_pair(2000004, -2000005)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1330, PID_pair(-2000004, 2000005)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1331, PID_pair(2000005, -2000005)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1332, PID_pair(2000006, -2000006)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1333, PID_pair(2000006, -1000001)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1333, PID_pair(-2000006, 1000001)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1334, PID_pair(2000006, -1000003)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1334, PID_pair(-2000006, 1000003)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1335, PID_pair(2000006, -1000005)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1335, PID_pair(-2000006, 1000005)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1338, PID_pair(2000006, -2000005)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1338, PID_pair(-2000006, 2000005)) );
       
        // Squark--squark production
        process_to_PIDs.insert( std::pair<int,PID_pair>(1351, PID_pair(1000001, 1000001)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1351, PID_pair(-1000001, -1000001)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1352, PID_pair(1000001, 1000003)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1352, PID_pair(-1000001, -1000003)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1353, PID_pair(1000001, 1000005)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1353, PID_pair(-1000001, -1000005)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1354, PID_pair(1000001, 2000001)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1354, PID_pair(-1000001, -2000001)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1355, PID_pair(1000001, 2000003)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1355, PID_pair(-1000001, -2000003)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1356, PID_pair(1000001, 2000005)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1356, PID_pair(-1000001, -2000005)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1357, PID_pair(1000002, 1000002)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1357, PID_pair(-1000002, -1000002)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1358, PID_pair(1000002, 1000004)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1358, PID_pair(-1000002, -1000004)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1360, PID_pair(1000002, 2000002)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1360, PID_pair(-1000002, -2000002)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1361, PID_pair(1000002, 2000004)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1361, PID_pair(-1000002, -2000004)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1363, PID_pair(1000002, 1000001)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1363, PID_pair(-1000002, -1000001)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1364, PID_pair(1000002, 1000003)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1364, PID_pair(-1000002, -1000003)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1365, PID_pair(1000002, 1000005)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1365, PID_pair(-1000002, -1000005)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1366, PID_pair(1000002, 2000001)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1366, PID_pair(-1000002, -2000001)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1367, PID_pair(1000002, 2000003)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1367, PID_pair(-1000002, -2000003)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1368, PID_pair(1000002, 2000005)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1368, PID_pair(-1000002, -2000005)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1369, PID_pair(1000003, 1000003)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1369, PID_pair(-1000003, -1000003)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1370, PID_pair(1000003, 1000005)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1370, PID_pair(-1000003, -1000005)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1371, PID_pair(1000003, 2000001)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1371, PID_pair(-1000003, -2000001)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1372, PID_pair(1000003, 2000003)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1372, PID_pair(-1000003, -2000003)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1373, PID_pair(1000003, 2000005)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1373, PID_pair(-1000003, -2000005)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1374, PID_pair(1000004, 1000004)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1374, PID_pair(-1000004, -1000004)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1376, PID_pair(1000004, 2000002)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1376, PID_pair(-1000004, -2000002)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1377, PID_pair(1000004, 2000004)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1377, PID_pair(-1000004, -2000004)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1379, PID_pair(1000004, 1000001)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1379, PID_pair(-1000004, -1000001)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1380, PID_pair(1000004, 1000003)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1380, PID_pair(-1000004, -1000003)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1381, PID_pair(1000004, 1000005)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1381, PID_pair(-1000004, -1000005)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1382, PID_pair(1000004, 2000001)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1382, PID_pair(-1000004, -2000001)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1383, PID_pair(1000004, 2000003)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1383, PID_pair(-1000004, -2000003)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1384, PID_pair(1000004, 2000005)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1384, PID_pair(-1000004, -2000005)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1385, PID_pair(1000005, 1000005)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1385, PID_pair(-1000005, -1000005)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1386, PID_pair(1000005, 2000001)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1386, PID_pair(-1000005, -2000001)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1387, PID_pair(1000005, 2000003)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1387, PID_pair(-1000005, -2000003)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1388, PID_pair(1000005, 2000005)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1388, PID_pair(-1000005, -2000005)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1393, PID_pair(1000006, 1000001)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1393, PID_pair(-1000006, -1000001)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1394, PID_pair(1000006, 1000003)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1394, PID_pair(-1000006, -1000003)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1395, PID_pair(1000006, 1000005)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1395, PID_pair(-1000006, -1000005)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1398, PID_pair(1000006, 2000005)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1398, PID_pair(-1000006, -2000005)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1399, PID_pair(2000001, 2000001)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1399, PID_pair(-2000001, -2000001)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1400, PID_pair(2000001, 2000003)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1400, PID_pair(-2000001, -2000003)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1401, PID_pair(2000001, 2000005)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1401, PID_pair(-2000001, -2000005)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1402, PID_pair(2000002, 2000002)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1402, PID_pair(-2000002, -2000002)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1403, PID_pair(2000002, 2000004)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1403, PID_pair(-2000002, -2000004)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1405, PID_pair(2000002, 1000001)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1405, PID_pair(-2000002, -1000001)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1406, PID_pair(2000002, 1000003)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1406, PID_pair(-2000002, -1000003)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1407, PID_pair(2000002, 1000005)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1407, PID_pair(-2000002, -1000005)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1408, PID_pair(2000002, 2000001)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1408, PID_pair(-2000002, -2000001)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1409, PID_pair(2000002, 2000003)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1409, PID_pair(-2000002, -2000003)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1410, PID_pair(2000002, 2000005)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1410, PID_pair(-2000002, -2000005)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1411, PID_pair(2000003, 2000003)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1411, PID_pair(-2000003, -2000003)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1412, PID_pair(2000003, 2000005)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1412, PID_pair(-2000003, -2000005)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1413, PID_pair(2000004, 2000004)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1413, PID_pair(-2000004, -2000004)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1415, PID_pair(2000004, 1000001)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1415, PID_pair(-2000004, -1000001)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1416, PID_pair(2000004, 1000003)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1416, PID_pair(-2000004, -1000003)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1417, PID_pair(2000004, 1000005)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1417, PID_pair(-2000004, -1000005)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1418, PID_pair(2000004, 2000001)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1418, PID_pair(-2000004, -2000001)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1419, PID_pair(2000004, 2000003)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1419, PID_pair(-2000004, -2000003)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1420, PID_pair(2000004, 2000005)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1420, PID_pair(-2000004, -2000005)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1421, PID_pair(2000005, 2000005)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1421, PID_pair(-2000005, -2000005)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1423, PID_pair(2000006, 1000001)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1423, PID_pair(-2000006, -1000001)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1424, PID_pair(2000006, 1000003)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1424, PID_pair(-2000006, -1000003)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1425, PID_pair(2000006, 1000005)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1425, PID_pair(-2000006, -1000005)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1428, PID_pair(2000006, 2000005)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1428, PID_pair(-2000006, -2000005)) );
        
        // Squark--electroweakino associated production
        process_to_PIDs.insert( std::pair<int,PID_pair>(1431, PID_pair(1000022, 1000001)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1431, PID_pair(1000022, -1000001)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1432, PID_pair(1000022, 1000002)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1432, PID_pair(1000022, -1000002)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1433, PID_pair(1000022, 1000003)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1433, PID_pair(1000022, -1000003)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1434, PID_pair(1000022, 1000004)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1434, PID_pair(1000022, -1000004)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1435, PID_pair(1000022, 1000005)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1435, PID_pair(1000022, -1000005)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1437, PID_pair(1000022, 2000001)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1437, PID_pair(1000022, -2000001)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1438, PID_pair(1000022, 2000002)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1438, PID_pair(1000022, -2000002)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1439, PID_pair(1000022, 2000003)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1439, PID_pair(1000022, -2000003)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1440, PID_pair(1000022, 2000004)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1440, PID_pair(1000022, -2000004)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1441, PID_pair(1000022, 2000005)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1441, PID_pair(1000022, -2000005)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1443, PID_pair(1000023, 1000001)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1443, PID_pair(1000023, -1000001)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1444, PID_pair(1000023, 1000002)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1444, PID_pair(1000023, -1000002)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1445, PID_pair(1000023, 1000003)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1445, PID_pair(1000023, -1000003)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1446, PID_pair(1000023, 1000004)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1446, PID_pair(1000023, -1000004)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1447, PID_pair(1000023, 1000005)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1447, PID_pair(1000023, -1000005)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1449, PID_pair(1000023, 2000001)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1449, PID_pair(1000023, -2000001)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1450, PID_pair(1000023, 2000002)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1450, PID_pair(1000023, -2000002)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1451, PID_pair(1000023, 2000003)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1451, PID_pair(1000023, -2000003)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1452, PID_pair(1000023, 2000004)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1452, PID_pair(1000023, -2000004)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1453, PID_pair(1000023, 2000005)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1453, PID_pair(1000023, -2000005)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1455, PID_pair(1000025, 1000001)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1455, PID_pair(1000025, -1000001)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1456, PID_pair(1000025, 1000002)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1456, PID_pair(1000025, -1000002)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1457, PID_pair(1000025, 1000003)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1457, PID_pair(1000025, -1000003)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1458, PID_pair(1000025, 1000004)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1458, PID_pair(1000025, -1000004)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1459, PID_pair(1000025, 1000005)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1459, PID_pair(1000025, -1000005)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1461, PID_pair(1000025, 2000001)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1461, PID_pair(1000025, -2000001)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1462, PID_pair(1000025, 2000002)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1462, PID_pair(1000025, -2000002)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1463, PID_pair(1000025, 2000003)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1463, PID_pair(1000025, -2000003)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1464, PID_pair(1000025, 2000004)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1464, PID_pair(1000025, -2000004)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1465, PID_pair(1000025, 2000005)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1465, PID_pair(1000025, -2000005)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1467, PID_pair(1000035, 1000001)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1467, PID_pair(1000035, -1000001)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1468, PID_pair(1000035, 1000002)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1468, PID_pair(1000035, -1000002)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1469, PID_pair(1000035, 1000003)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1469, PID_pair(1000035, -1000003)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1470, PID_pair(1000035, 1000004)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1470, PID_pair(1000035, -1000004)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1471, PID_pair(1000035, 1000005)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1471, PID_pair(1000035, -1000005)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1473, PID_pair(1000035, 2000001)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1473, PID_pair(1000035, -2000001)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1474, PID_pair(1000035, 2000002)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1474, PID_pair(1000035, -2000002)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1475, PID_pair(1000035, 2000003)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1475, PID_pair(1000035, -2000003)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1476, PID_pair(1000035, 2000004)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1476, PID_pair(1000035, -2000004)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1477, PID_pair(1000035, 2000005)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1477, PID_pair(1000035, -2000005)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1491, PID_pair(1000024, 1000001)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1491, PID_pair(-1000024, -1000001)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1492, PID_pair(-1000024, 1000002)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1492, PID_pair(1000024, -1000002)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1493, PID_pair(1000024, 1000003)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1493, PID_pair(-1000024, -1000003)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1494, PID_pair(-1000024, 1000004)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1494, PID_pair(1000024, -1000004)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1495, PID_pair(1000024, 1000005)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1495, PID_pair(-1000024, -1000005)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1496, PID_pair(-1000024, 1000006)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1496, PID_pair(1000024, -1000006)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1501, PID_pair(1000024, 2000005)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1501, PID_pair(-1000024, -2000005)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1502, PID_pair(-1000024, 2000006)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1502, PID_pair(1000024, -2000006)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1503, PID_pair(1000037, 1000001)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1503, PID_pair(-1000037, -1000001)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1504, PID_pair(-1000037, 1000002)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1504, PID_pair(1000037, -1000002)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1505, PID_pair(1000037, 1000003)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1505, PID_pair(-1000037, -1000003)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1506, PID_pair(-1000037, 1000004)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1506, PID_pair(1000037, -1000004)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1507, PID_pair(1000037, 1000005)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1507, PID_pair(-1000037, -1000005)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1508, PID_pair(-1000037, 1000006)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1508, PID_pair(1000037, -1000006)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1513, PID_pair(1000037, 2000005)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1513, PID_pair(-1000037, -2000005)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1514, PID_pair(-1000037, 2000006)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1514, PID_pair(1000037, -2000006)) );
        
        // Neutralino & chargino production
        process_to_PIDs.insert( std::pair<int,PID_pair>(1551, PID_pair(1000022, 1000022)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1552, PID_pair(1000022, 1000023)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1553, PID_pair(1000023, 1000023)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1554, PID_pair(1000022, 1000025)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1555, PID_pair(1000023, 1000025)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1556, PID_pair(1000025, 1000025)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1557, PID_pair(1000022, 1000035)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1558, PID_pair(1000023, 1000035)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1559, PID_pair(1000025, 1000035)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1560, PID_pair(1000035, 1000035)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1571, PID_pair(1000024, 1000022)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1572, PID_pair(-1000024, 1000022)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1573, PID_pair(1000037, 1000022)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1574, PID_pair(-1000037, 1000022)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1575, PID_pair(1000024, 1000023)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1576, PID_pair(-1000024, 1000023)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1577, PID_pair(1000037, 1000023)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1578, PID_pair(-1000037, 1000023)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1579, PID_pair(1000024, 1000025)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1580, PID_pair(-1000024, 1000025)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1581, PID_pair(1000037, 1000025)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1582, PID_pair(-1000037, 1000025)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1583, PID_pair(1000024, 1000035)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1584, PID_pair(-1000024, 1000035)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1585, PID_pair(1000037, 1000035)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1586, PID_pair(-1000037, 1000035)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1591, PID_pair(1000024, -1000024)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1592, PID_pair(1000024, -1000037)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1593, PID_pair(1000037, -1000024)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1594, PID_pair(1000037, -1000037)) );
        
        // Gluino--electroweakino associated production
        process_to_PIDs.insert( std::pair<int,PID_pair>(1601, PID_pair(1000021, 1000022)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1602, PID_pair(1000021, 1000023)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1603, PID_pair(1000021, 1000025)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1604, PID_pair(1000021, 1000035)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1621, PID_pair(1000021, 1000024)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1621, PID_pair(1000021, -1000024)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1622, PID_pair(1000021, 1000037)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1622, PID_pair(1000021, -1000037)) );
        
        // Slepton production
        process_to_PIDs.insert( std::pair<int,PID_pair>(1651, PID_pair(1000011, -1000011)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1662, PID_pair(1000012, -1000012)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1673, PID_pair(1000012, -1000011)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1673, PID_pair(-1000012, 1000011)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1679, PID_pair(1000013, -1000013)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1688, PID_pair(1000014, -1000014)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1698, PID_pair(1000014, -1000013)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1698, PID_pair(-1000014, 1000013)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1703, PID_pair(1000015, -1000015)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1708, PID_pair(1000015, -2000015)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1709, PID_pair(2000015, -1000015)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1710, PID_pair(1000016, -1000016)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1719, PID_pair(1000016, -1000015)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1719, PID_pair(-1000016, 1000015)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1722, PID_pair(1000016, -2000015)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1722, PID_pair(-1000016, 2000015)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1723, PID_pair(2000011, -2000011)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1739, PID_pair(2000013, -2000013)) );
        process_to_PIDs.insert( std::pair<int,PID_pair>(1751, PID_pair(2000015, -2000015)) );

        first = false;
      }


      if (*Loop::iteration == COLLIDER_INIT)
      {
        result.clear();
      }

      if (*Loop::iteration == XSEC_CALCULATION)
      {
        std::vector<int> process_codes = *Dep::ProcessCodes;
        for (int& c : process_codes)
        {
          // TODO: Should we use multimap::find here, which only finds the first
          // matching process code, or should we use multimap::equal_range, 
          // to get all the matching elements? 
          // Depends on what the rest of the code assumes...
          auto it = process_to_PIDs.find(c);
          if (it == process_to_PIDs.end())
          {
            std::stringstream errmsg_ss;
            errmsg_ss << "Can't find the Pythia process code " << c << " in the process_to_PIDs map." << endl;
            ColliderBit_error().raise(LOCAL_INFO, errmsg_ss.str());
          } 
          else
          {
            PID_pair& p = process_to_PIDs.find(c)->second;
            result.push_back(p);
          }
        }
      }

      // _Anders
      if (*Loop::iteration == START_SUBPROCESS)
      {
        cout << DEBUG_PREFIX << "getProcessPIDPairs: it = START_SUBPROCESS, result.size() = " << result.size() << endl;
        for (PID_pair& p : result)
        {
          cout << DEBUG_PREFIX << "getProcessPIDPairs: it = START_SUBPROCESS, result element: (" << p.first << "," << p.second << ")" << endl;          
        }
      }

    }

  } 
} 


