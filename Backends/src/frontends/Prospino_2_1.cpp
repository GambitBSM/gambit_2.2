//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  Frontend for Prospino 2.1 backend
///
///  *********************************************
///
///  Authors (add name and date if you modify):
///
///  \author Anders Kvellestad
///  \date 2019 Apr
///
///  *********************************************

#include "gambit/Backends/backend_singleton.hpp"
#include "gambit/Backends/frontend_macros.hpp"
#include "gambit/Backends/frontends/Prospino_2_1.hpp"
#include "gambit/Elements/mssm_slhahelp.hpp"
#include "gambit/Elements/slhaea_helpers.hpp"

#include "gambit/Utils/version.hpp"

#define BACKEND_DEBUG 0


// Backend init function
BE_INI_FUNCTION
{
    // Help Prospino find itself
    std::string prospino_dir = Backends::backendInfo().path_dir(STRINGIFY(BACKENDNAME), STRINGIFY(VERSION));
    Fstring<500> prospino_dir_in = prospino_dir.c_str();
    prospino_gb_init(prospino_dir_in);
}
END_BE_INI_FUNCTION


// Convenience functions (definition)
BE_NAMESPACE
{
  // Convenience function to run Prospino and get a vector of cross-sections
  map_str_dbl run_prospino(const SLHAstruct& slha, prospino_settings& ps)
  {

    // Get type converter 
    using SLHAea::to;

    std::cout << "DEBUG: run_prospino: Begin..." << std::endl;

    Farray<Fdouble,1,20> unimass;
    Farray<Fdouble,0,99> lowmass;

    Farray<Fdouble,1,2,1,2> uu_in, vv_in;
    Farray<Fdouble,1,4,1,4> bw_in;
    Farray<Fdouble,1,2,1,2> mst_in, msb_in, msl_in;


    lowmass(0) = to<double>(slha.at("HMIX").at(1).at(1));
    lowmass(1) = to<double>(slha.at("MSOFT").at(1).at(1));
    lowmass(2) = to<double>(slha.at("MSOFT").at(2).at(1));
    lowmass(3) = to<double>(slha.at("MSOFT").at(3).at(1));

    lowmass(4) = to<double>(slha.at("MASS").at(1000021).at(1));

    lowmass(5) = to<double>(slha.at("MASS").at(1000022).at(1));
    lowmass(6) = to<double>(slha.at("MASS").at(1000023).at(1));
    lowmass(7) = to<double>(slha.at("MASS").at(1000025).at(1));
    lowmass(8) = to<double>(slha.at("MASS").at(1000035).at(1));

    lowmass(9) = to<double>(slha.at("MASS").at(1000024).at(1));
    lowmass(10) = to<double>(slha.at("MASS").at(1000037).at(1));

    lowmass(11) = to<double>(slha.at("MASS").at(1000001).at(1));
    lowmass(52) = lowmass(11);

    lowmass(12) = to<double>(slha.at("MASS").at(2000001).at(1));
    lowmass(58) = lowmass(12);

    lowmass(13) = to<double>(slha.at("MASS").at(1000002).at(1));
    lowmass(51) = lowmass(13);

    lowmass(14) = to<double>(slha.at("MASS").at(2000002).at(1));
    lowmass(57) = lowmass(14);

    lowmass(17) = to<double>(slha.at("MASS").at(1000005).at(1));
    lowmass(55) = lowmass(17);

    lowmass(18) = to<double>(slha.at("MASS").at(2000005).at(1));
    lowmass(61) = lowmass(18);

    lowmass(19) = to<double>(slha.at("MASS").at(1000006).at(1));
    lowmass(56) = lowmass(19);

    lowmass(20) = to<double>(slha.at("MASS").at(2000006).at(1));
    lowmass(62) = lowmass(20);

    lowmass(30) = to<double>(slha.at("MASS").at(1000011).at(1));

    lowmass(31) = to<double>(slha.at("MASS").at(2000011).at(1));

    lowmass(32) = to<double>(slha.at("MASS").at(1000012).at(1));

    lowmass(33) = to<double>(slha.at("MASS").at(1000015).at(1));

    lowmass(34) = to<double>(slha.at("MASS").at(2000015).at(1));

    lowmass(35) = to<double>(slha.at("MASS").at(1000016).at(1));

    lowmass(40) = to<double>(slha.at("MASS").at(36).at(1));

    lowmass(41) = to<double>(slha.at("MASS").at(25).at(1));

    lowmass(42) = to<double>(slha.at("MASS").at(35).at(1));

    lowmass(43) = to<double>(slha.at("MASS").at(37).at(1));

    lowmass(53) = to<double>(slha.at("MASS").at(1000003).at(1));

    lowmass(59) = to<double>(slha.at("MASS").at(2000003).at(1));

    lowmass(54) = to<double>(slha.at("MASS").at(1000004).at(1));

    lowmass(60) = to<double>(slha.at("MASS").at(2000004).at(1));

    // Degenerate squark masses
    double degen_squark_mass_8 = 0.;
    degen_squark_mass_8 += lowmass(51) + lowmass(52) + lowmass(53) + lowmass(54);
    degen_squark_mass_8 += lowmass(57) + lowmass(58) + lowmass(59) + lowmass(60);
    degen_squark_mass_8 = degen_squark_mass_8/8.;
    lowmass(15) = degen_squark_mass_8;

    double degen_squark_mass_10 = 0.;
    degen_squark_mass_10 += lowmass(51) + lowmass(52) + lowmass(53) + lowmass(54) + lowmass(55);
    degen_squark_mass_10 += lowmass(57) + lowmass(58) + lowmass(59) + lowmass(60) + lowmass(61);
    degen_squark_mass_10 = degen_squark_mass_10/10.;
    lowmass(16) = degen_squark_mass_10;

    // BLOCK ALPHA
    double alpha = to<double>( slha.at("ALPHA").back().front() );
    lowmass(44) = sin(alpha);
    lowmass(45) = cos(alpha);

    // AD, AU, AE
    lowmass(21) = -to<double>(slha.at("AD").at(3,3).at(2));   // Note sign!
    lowmass(24) = -to<double>(slha.at("AU").at(3,3).at(2));   // Note sign!
    lowmass(36) = -to<double>(slha.at("AE").at(3,3).at(2));   // Note sign!

    // Neutralino mixing matrix
    bw_in(1,1) = to<double>(slha.at("NMIX").at(1,1).at(2));
    bw_in(1,2) = to<double>(slha.at("NMIX").at(1,2).at(2));
    bw_in(1,3) = to<double>(slha.at("NMIX").at(1,3).at(2));
    bw_in(1,4) = to<double>(slha.at("NMIX").at(1,4).at(2));
    bw_in(2,1) = to<double>(slha.at("NMIX").at(2,1).at(2));
    bw_in(2,2) = to<double>(slha.at("NMIX").at(2,2).at(2));
    bw_in(2,3) = to<double>(slha.at("NMIX").at(2,3).at(2));
    bw_in(2,4) = to<double>(slha.at("NMIX").at(2,4).at(2));
    bw_in(3,1) = to<double>(slha.at("NMIX").at(3,1).at(2));
    bw_in(3,2) = to<double>(slha.at("NMIX").at(3,2).at(2));
    bw_in(3,3) = to<double>(slha.at("NMIX").at(3,3).at(2));
    bw_in(3,4) = to<double>(slha.at("NMIX").at(3,4).at(2));
    bw_in(4,1) = to<double>(slha.at("NMIX").at(4,1).at(2));
    bw_in(4,2) = to<double>(slha.at("NMIX").at(4,2).at(2));
    bw_in(4,3) = to<double>(slha.at("NMIX").at(4,3).at(2));
    bw_in(4,4) = to<double>(slha.at("NMIX").at(4,4).at(2));

    // Chargino mixing matrix
    uu_in(1,1) = to<double>(slha.at("UMIX").at(1,1).at(2));
    uu_in(1,2) = to<double>(slha.at("UMIX").at(1,2).at(2));
    uu_in(2,1) = to<double>(slha.at("UMIX").at(2,1).at(2));
    uu_in(2,2) = to<double>(slha.at("UMIX").at(2,2).at(2));

    vv_in(1,1) = to<double>(slha.at("VMIX").at(1,1).at(2));
    vv_in(1,2) = to<double>(slha.at("VMIX").at(1,2).at(2));
    vv_in(2,1) = to<double>(slha.at("VMIX").at(2,1).at(2));
    vv_in(2,2) = to<double>(slha.at("VMIX").at(2,2).at(2));

    // Sfermion mixing matrices
    mst_in(1,1) = to<double>(slha.at("STOPMIX").at(1,1).at(2));
    mst_in(1,2) = to<double>(slha.at("STOPMIX").at(1,2).at(2));
    mst_in(2,1) = to<double>(slha.at("STOPMIX").at(2,1).at(2));
    mst_in(2,2) = to<double>(slha.at("STOPMIX").at(2,2).at(2));

    msb_in(1,1) = to<double>(slha.at("SBOTMIX").at(1,1).at(2));
    msb_in(1,2) = to<double>(slha.at("SBOTMIX").at(1,2).at(2));
    msb_in(2,1) = to<double>(slha.at("SBOTMIX").at(2,1).at(2));
    msb_in(2,2) = to<double>(slha.at("SBOTMIX").at(2,2).at(2));

    msl_in(1,1) = to<double>(slha.at("STAUMIX").at(1,1).at(2));
    msl_in(1,2) = to<double>(slha.at("STAUMIX").at(1,2).at(2));
    msl_in(2,1) = to<double>(slha.at("STAUMIX").at(2,1).at(2));
    msl_in(2,2) = to<double>(slha.at("STAUMIX").at(2,2).at(2));


    // Call prospino
    Farray<Fdouble,0,6> prospino_result;

    prospino_gb(prospino_result, ps.inlo, ps.isq_ng_in, ps.icoll_in, ps.energy_in, ps.i_error_in, 
                ps.final_state_in, ps.ipart1_in, ps.ipart2_in, ps.isquark1_in, ps.isquark2_in,
                unimass, lowmass, uu_in, vv_in, bw_in, mst_in, msb_in, msl_in);

    
    std::cout << "DEBUG: run_prospino: " << std::endl;
    std::cout << "DEBUG: run_prospino: prospino_result(0) = " << prospino_result(0) << std::endl;
    std::cout << "DEBUG: run_prospino: prospino_result(1) = " << prospino_result(1) << std::endl;
    std::cout << "DEBUG: run_prospino: prospino_result(2) = " << prospino_result(2) << std::endl;
    std::cout << "DEBUG: run_prospino: prospino_result(3) = " << prospino_result(3) << std::endl;
    std::cout << "DEBUG: run_prospino: prospino_result(4) = " << prospino_result(4) << std::endl;
    std::cout << "DEBUG: run_prospino: prospino_result(5) = " << prospino_result(5) << std::endl;
    std::cout << "DEBUG: run_prospino: prospino_result(6) = " << prospino_result(6) << std::endl;
    std::cout << "DEBUG: run_prospino: " << std::endl;
    std::cout << "DEBUG: run_prospino: ...End" << std::endl;


    // Fill the result map with the content of prospino_result
    map_str_dbl result;
    result["LO[pb]"] = prospino_result(0);
    result["LO_rel_error"] = prospino_result(1);
    result["NLO[pb]"] = prospino_result(2);
    result["NLO_rel_error"] = prospino_result(3);
    result["K"] = prospino_result(4);
    result["LO_ms[pb]"] = prospino_result(5);
    result["NLO_ms[pb]"] = prospino_result(6);

    return result;

  }
}
END_BE_NAMESPACE
