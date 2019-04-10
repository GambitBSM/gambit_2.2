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

#include "gambit/Elements/spectrum.hpp"
// #include "gambit/Elements/spectrum_factories.hpp"
// #include "gambit/Models/SimpleSpectra/NMSSMSimpleSpec.hpp"

#include "gambit/Utils/version.hpp"

#define BACKEND_DEBUG 0


// Convenience functions (definition)
BE_NAMESPACE
{
  // Convenience function to run SPheno and obtain the spectrum
  std::vector<double> run_prospino(const Spectrum& spectrum)
  {

    std::cout << "DEBUG: run_prospino: Begin..." << std::endl;

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

// Initialisation function (definition)
BE_INI_FUNCTION{} END_BE_INI_FUNCTION
