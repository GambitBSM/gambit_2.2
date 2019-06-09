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
  map_str_dbl run_prospino(const SLHAstruct& slha_in, const param_map_type& params, prospino_settings& ps)
  {

    // Get type converter 
    using SLHAea::to;

    std::cout << "DEBUG: run_prospino: Begin..." << std::endl;

    // Copy the slha object so that we can modify it
    SLHAstruct slha(slha_in);

    // Contstruct EXTPAR block from input parameters
    SLHAea_add_block(slha, "EXTPAR");

    slha["EXTPAR"][""] << 0 << *params.at("Qin") << "# scale Q where the parameters below are defined";

    slha["EXTPAR"][""] << 1 << *params.at("M1") << "# M_1";
    slha["EXTPAR"][""] << 2 << *params.at("M2") << "# M_2";
    slha["EXTPAR"][""] << 3 << *params.at("M3") << "# M_3";

    slha["EXTPAR"][""] << 11 << *params.at("Au_33") << "# A_t";
    slha["EXTPAR"][""] << 12 << *params.at("Ad_33") << "# A_b";
    slha["EXTPAR"][""] << 13 << *params.at("Ae_33") << "# A_l";

    // slha["EXTPAR"][""] << 21 << *params.at("mHd2") << "# m_Hd^2";
    // slha["EXTPAR"][""] << 22 << *params.at("mHd2") << "# m_Hu^2";
    slha["EXTPAR"][""] << 23 << *params.at("mu") << "# mu";
    slha["EXTPAR"][""] << 24 << pow(*params.at("mA"),2) << "# m_A^2";

    slha["EXTPAR"][""] << 31 << sqrt(*params.at("ml2_11")) << "# M_(L,11)";
    slha["EXTPAR"][""] << 32 << sqrt(*params.at("ml2_22")) << "# M_(L,22)";
    slha["EXTPAR"][""] << 33 << sqrt(*params.at("ml2_33")) << "# M_(L,33)";
    slha["EXTPAR"][""] << 34 << sqrt(*params.at("me2_11")) << "# M_(E,11)";
    slha["EXTPAR"][""] << 35 << sqrt(*params.at("me2_22")) << "# M_(E,22)";
    slha["EXTPAR"][""] << 36 << sqrt(*params.at("me2_33")) << "# M_(E,33)";
    slha["EXTPAR"][""] << 41 << sqrt(*params.at("mq2_11")) << "# M_(Q,11)";
    slha["EXTPAR"][""] << 42 << sqrt(*params.at("mq2_22")) << "# M_(Q,22)";
    slha["EXTPAR"][""] << 43 << sqrt(*params.at("mq2_33")) << "# M_(Q,33)";
    slha["EXTPAR"][""] << 44 << sqrt(*params.at("mu2_11")) << "# M_(U,11)";
    slha["EXTPAR"][""] << 45 << sqrt(*params.at("mu2_22")) << "# M_(U,22)";
    slha["EXTPAR"][""] << 46 << sqrt(*params.at("mu2_33")) << "# M_(U,33)";
    slha["EXTPAR"][""] << 47 << sqrt(*params.at("md2_11")) << "# M_(D,11)";
    slha["EXTPAR"][""] << 48 << sqrt(*params.at("md2_22")) << "# M_(D,22)";
    slha["EXTPAR"][""] << 49 << sqrt(*params.at("md2_33")) << "# M_(D,33)";

    std::cout << "DEBUG: SLHAstruct content:" << std::endl;
    std::cout << slha.str() << std::endl;
    

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
