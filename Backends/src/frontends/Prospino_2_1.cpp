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

#include "gambit/Backends/frontend_macros.hpp"
#include "gambit/Backends/frontends/Prospino_2_1.hpp"
#include "gambit/Elements/mssm_slhahelp.hpp"
#include "gambit/Elements/slhaea_helpers.hpp"

#include <typeinfo>

// #include "gambit/Elements/spectrum.hpp"
// #include "gambit/Elements/spectrum_factories.hpp"
// #include "gambit/Models/SimpleSpectra/NMSSMSimpleSpec.hpp"

#include "gambit/Utils/version.hpp"

#define BACKEND_DEBUG 0


// Convenience functions (definition)
BE_NAMESPACE
{
  // Convenience function to run Prospino and get a vector of cross-sections
  // std::vector<double> run_prospino(const Spectrum& spectrum)
  std::vector<double> run_prospino(const SLHAstruct& slha_in, const param_map_type& params)
  {

    // // Get path
    // // TODO: move this to init function
    // const str be = "Pythia" + model_suffix;
    // const str ver = Backends::backendInfo().default_version(be);
    // pythia_doc_path = Backends::backendInfo().path_dir(be, ver) + "/../share/Pythia8/xmldoc/";


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

    // std::cout << "DEBUG:  slha.at('MASS').at(25).at(1) = " << to<double>(slha.at("MASS").at(25).at(1)) << std::endl;
    // std::cout << "DEBUG:  slha.at('EXTPAR').at(1).at(1) = " << to<double>(slha.at("EXTPAR").at(1).at(1)) << std::endl;
    
    std::cout << "DEBUG: SLHAstruct content:" << std::endl;
    std::cout << slha.str() << std::endl;
    


    // std::cout << "DEBUG:  slha.at('EXTPAR').at(1).at(1) = " << to<double>(slha.at("EXTPAR").at(1).at(1)) << std::endl;
    // std::cout << "DEBUG:  slha.at('MASS').at(25).at(1) = " << to<double>(slha.at("MASS").at(25).at(1)) << std::endl;


    Finteger inlo = 1;                 // specify LO only[0] or complete NLO (slower)[1]
    Finteger isq_ng_in = 1;            // specify degenerate [0] or free [1] squark masses
    Finteger icoll_in = 1;             // collider : tevatron[0], lhc[1]
    Fdouble energy_in = 13000.0;       // collider energy in GeV
    Finteger i_error_in = 0;           // with central scale [0] or scale variation [1]

    Fstring<2> final_state_in = "nn";  // select process
    Finteger ipart1_in = 1;            //
    Finteger ipart2_in = 2;            //
    Finteger isquark1_in = 0;          //
    Finteger isquark2_in = 0;          //

    Farray<Fdouble,1,20> unimass;
    Farray<Fdouble,0,99> lowmass;

    Farray<Fdouble,1,2,1,2> uu_in, vv_in;
    Farray<Fdouble,1,4,1,4> bw_in;
    Farray<Fdouble,1,2,1,2> mst_in, msb_in, msl_in;

/*
  integer                            :: ipart1, ipart2   ! set in INIT_SUSY   : internal variable for final state particles 
  real(kind=double)                  :: mg,ms            ! set in INIT_SUSY
  real(kind=double)                  :: mh1,mh2,mch      ! set in INIT_SUSY
  real(kind=double)                  :: sin_a,cos_a      ! set in INIT_SUSY
  real(kind=double)                  :: mu_susy,tan_b    ! set in INIT_SUSY
  real(kind=double)                  :: a_b,a_t          ! set in INIT_SUSY
  real(kind=double), dimension(-6:6) :: msq              ! set in INIT_SUSY
  real(kind=double), dimension(1:8)  :: smass_n,mass_n   ! set in INIT_SUSY
  real(kind=double), dimension(2,2)  :: uu,vv            ! set in INIT_SUSY
  real(kind=double), dimension(2,2)  :: mst,msb,msl      ! set in INIT_SUSY
  real(kind=double), dimension(4,4)  :: bw               ! set in INIT_SUSY
  real(kind=double), dimension(1:4)  :: mass_s           ! set in INIT_SUSY
  real(kind=double), dimension(1:4)  :: mass_x           ! set in INIT_SUSY
  real(kind=double)                  :: mg_orig,ms_orig  ! set in INIT_SUSY
  complex(kind=double), dimension(4,4) :: zz             ! set in INIT_SUSY
*/


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

/*
  output parameters: 
  lowmass(0)  mu                                 c
  lowmass(1)  m_1
  lowmass(2)  m_2
  lowmass(3)  m_3
                                                 
  lowmass(4)  gluino mass
  lowmass(5)  \
  lowmass(6)   \
  lowmass(7)   /  neutralino masses [with sign]
  lowmass(8)  /
  lowmass(9)  \
  lowmass(10) /   chargino masses
                                                 
  lowmass(15) degenerate squark mass (8)
  lowmass(16) degenerate squark mass (10)
                                                 
  lowmass(21) a_b
  lowmass(24) a_t
                                                 
  lowmass(30) selectron_l mass
  lowmass(31) selectron_r mass
  lowmass(32) selectron-neutrino mass
  lowmass(33) stau_1 mass
  lowmass(34) stau_2 mass
  lowmass(35) stau-neutrino mass
  lowmass(36) a_tau
                                                 
  lowmass(40) cp odd higgs mass
  lowmass(41) light cp even higgs mass
  lowmass(42) heavy cp even higgs mass
  lowmass(43) charged higgs mass
  lowmass(44) sin(alpha
  lowmass(45) cos(alpha
                                                 
  like cteq: u,d,s,c,b,t first L then R
  lowmass(51) sup_L
  lowmass(52) sdown_L
  lowmass(53) sstrange_L
  lowmass(54) scharm_L
  lowmass(55) sbottom_1
  lowmass(56) stop_1
  lowmass(57) sup_R
  lowmass(58) sdown_R
  lowmass(59) sstrange_R
  lowmass(60) scharm_R
  lowmass(61) sbottom_2
  lowmass(62) stop_2
                                                 
  lowmass(80) unification scale
  lowmass(81) unified coupling alpha(m_x)
                                                 
  lowmass(91) trilinear higgs coupling lambda(1)
  .......
  lowmass(97) lambda(7)
                                                 
  lowmass(99
*/


    // (*MCha)(i) = spectrum.get(Par::Pole_Mass, "~chi+",i);
    // (*ZD)(i,j) = spectrum.get(Par::Pole_Mixing, "~d", i, j);


    // Call prospino
    prospino_gb(inlo, isq_ng_in, icoll_in, energy_in, i_error_in, final_state_in, ipart1_in, ipart2_in, isquark1_in, isquark2_in,
                unimass, lowmass, uu_in, vv_in, bw_in, mst_in, msb_in, msl_in);

    std::cout << "DEBUG: run_prospino: ...End" << std::endl;


    /*
      call PROSPINO_OPEN_CLOSE(0)                                                            ! open all input/output files
      
      call PROSPINO_CHECK_HIGGS(final_state_in)                                              ! lock Higgs final states
      call PROSPINO_CHECK_FS(final_state_in,ipart1_in,ipart2_in,lfinal)                      ! check final state 
      if (.not. lfinal ) then
         print*, " final state not correct ",final_state_in,ipart1_in,ipart2_in
         call HARD_STOP                                                                      ! finish if final state bad
      end if

      call PROSPINO(inlo,isq_ng_in,icoll_in,energy_in,i_error_in,final_state_in,ipart1_in,ipart2_in,isquark1_in,isquark2_in) ! actual prospino call
    */


    /* 
      (From the Prospino code documentation)

      Options for final_state_in:

      ng     neutralino/chargino + gluino  
      ns     neutralino/chargino + squark  
      nn     neutralino/chargino pair combinations  
      ll     slepton pair combinations  
      sb     squark-antisquark  
      ss     squark-squark  
      tb     stop-antistop  
      bb     sbottom-antisbottom  
      gg     gluino pair  
      sg     squark + gluino  
      lq     leptoquark pairs (using stop1 mass) 
      le     leptoquark plus lepton (using stop1 mass) 
      hh     charged Higgs pairs (private code only!)
      ht     charged Higgs with top (private code only!)

      Squark and antisquark added, but taking into account different sb or ss.


      Options for ipart1_in, ipart2_in:

      final_state_in = ng,ns,nn
      ipart1_in   = 1,2,3,4  neutralinos
                    5,6      positive charge charginos
                    7,8      negative charge charginos
      ipart2_in the same
          chargino+ and chargino- different processes
                                                                               
      final_state_in = ll
      ipart1_in   = 0        sel,sel + ser,ser  (first generation)
                    1        sel,sel
                    2        ser,ser
                    3        snel,snel
                    4        sel+,snl
                    5        sel-,snl
                    6        stau1,stau1
                    7        stau2,stau2
                    8        stau1,stau2
                    9        sntau,sntau
                   10        stau1+,sntau
                   11        stau1-,sntau
                   12        stau2+,sntau
                   13        stau2-,sntau
                   14        H+,H- in Drell-Yan channel
                                                                               
      final_state_in = tb and bb
      ipart1_in   = 1        stop1/sbottom1 pairs
                    2        stop2/sbottom2 pairs
                                                                               
      Note: Otherwise ipart1_in,ipart2_in have to set to one if not used.


      Options for isquark1_in, isquark1_in:

      for LO with light-squark flavor in the final state
      isquark1_in     =  -5,-4,-3,-2,-1,+1,+2,+3,+4,+5
                        (bL cL sL dL uL uR dR sR cR bR) in CteQ ordering
      isquark1_in     = 0 sum over light-flavor squarks throughout
                          (the squark mass in the data files is then averaged)

      flavors in initial state: only light-flavor partons, no bottoms
                                bottom partons only for Higgs channels

      flavors in final state: light-flavor quarks summed over five flavors
    */


    // Dummy result
    std::vector<double> xsec_vals;
    xsec_vals.push_back(1.2345);
    xsec_vals.push_back(1.2345 * 0.1);

    return xsec_vals;
  }
}
END_BE_NAMESPACE

// // Initialisation function (definition)
// BE_INI_FUNCTION{} END_BE_INI_FUNCTION

// Backend init function
BE_INI_FUNCTION
{
    Fstring<500> prospino_dir_in = "/home/anders/physics/GAMBIT/gambit/Backends/installed/prospino/2.1";
    prospino_gb_init(prospino_dir_in);
}
END_BE_INI_FUNCTION
