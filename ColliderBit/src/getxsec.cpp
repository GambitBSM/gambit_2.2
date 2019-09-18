//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  Module functions for computing cross-sections
///
///  *********************************************
///
///  Authors (add name and date if you modify):
///
///  \author Pat Scott
///          (p.scott@imperial.ac.uk)
///  \date 2019 Feb, May
///
///  \author Anders Kvellestad
///          (a.kvellestad@imperial.ac.uk)
///  \date 2019 Sep
///
///  *********************************************

#include "gambit/ColliderBit/ColliderBit_eventloop.hpp"
#include "gambit/ColliderBit/all_process_codes_to_PID_pairs.hpp"

#define COLLIDERBIT_DEBUG
#define DEBUG_PREFIX "DEBUG: OMP thread " << omp_get_thread_num() << ":  "

namespace Gambit
{

  namespace ColliderBit
  {

    // ======= Utility functions =======


    /// Helper function that takes a cross-section value in fb or pb, 
    /// along with an absolute or relative uncertainty, and returns the 
    /// xsec and absolute uncertainty in fb.
    std::pair<double,double> convert_xsecs_to_fb(double input_xsec, double input_xsec_uncert, str input_unit, bool input_fractional_uncert)
    {
      double xsec_fb;
      double xsec_uncert_fb;

      if (input_unit == "fb" && !input_fractional_uncert)
      {
        xsec_fb = input_xsec;
        xsec_uncert_fb = input_xsec_uncert;
      }
      else if (input_unit == "fb" && input_fractional_uncert)
      {
        xsec_fb = input_xsec;
        xsec_uncert_fb = input_xsec_uncert * xsec_fb;
      }
      else if (input_unit == "pb" && !input_fractional_uncert)
      {
        xsec_fb = input_xsec * 1000.;
        xsec_uncert_fb = input_xsec_uncert * 1000.;
      }
      else if (input_unit == "pb" && input_fractional_uncert)
      {
        xsec_fb = input_xsec * 1000.;
        xsec_uncert_fb = input_xsec_uncert * xsec_fb;
      }
      else
      {
        ColliderBit_error().raise(LOCAL_INFO, "Unknown combination of options for function convert_xsecs_to_fb.");
      }      

      return std::make_pair(xsec_fb, xsec_uncert_fb);
    }



    // ======= Module functions =======


    /// Dummy function for testing out how to return a PIDPairCrossSectionFunc
    xsec_container PIDPairCrossSection_dummy(const iipair& pids, const Spectrum& MSSM_spectrum)
    {
      xsec_container xs_result;

      // Get an SLHA1 object
      const SLHAstruct& slha = MSSM_spectrum.getSLHAea(1);

      // "Calculate" cross section
      // (only stop-stopbar production)
      double xs_fb = 0.0;
      if (pids.first == 1000002 && pids.second == -1000002)
      {
        xs_fb = 3.09816e-11 + 9.08223e-13;
      }
      else
      {
        xs_fb = 1.0;
      }
      double xs_rel_err = 0.01;
      // double xs_err = xs_fb * xs_rel_err;

      // Save result in xs_result
      xs_result.set_xsec(xs_fb, xs_rel_err);

      // Construct info string of the form "PID1:<PID1>, PID2:<PID2>"
      std::stringstream info_ss;
      info_ss << "PID1:" << pids.first << ", " << "PID2:" << pids.second;
      xs_result.set_info_string(info_ss.str());

      return xs_result;
    }





    /// Test functions for provding PIDPairCrossSectionsMap (cross-sections in fb)
    PID_pair_xsec_container silly_PID_xsec_constructor(double xsec_val)
    {
      PID_pair_xsec_container result;

      result.reset();
      result.set_xsec(xsec_val, xsec_val * 0.01);

      return result;
    }

    void getPIDPairCrossSectionsMap_testing(map_iipair_PID_pair_xsec& result)
    {
      using namespace Pipes::getPIDPairCrossSectionsMap_testing;

      static bool first = true;
      if (first)
      {
        // Gluino--gluino
        result[iipair(1000021,1000021)] = silly_PID_xsec_constructor(0.623E-01 * 1e3);
        // Gluino--squark
        result[iipair(1000001,1000021)] = silly_PID_xsec_constructor(0.14855 * 1e3);
        result[iipair(1000002,1000021)] = silly_PID_xsec_constructor(0.57564E-01 * 1e3);
        result[iipair(1000003,1000021)] = silly_PID_xsec_constructor(0.67117E-02 * 1e3);
        result[iipair(1000004,1000021)] = silly_PID_xsec_constructor(0.31759E-02 * 1e3);
        result[iipair(1000005,1000021)] = silly_PID_xsec_constructor(0 * 1e3);
        result[iipair(2000001,1000021)] = silly_PID_xsec_constructor(0.16430 * 1e3);
        result[iipair(2000002,1000021)] = silly_PID_xsec_constructor(0.65192E-01 * 1e3);
        result[iipair(2000003,1000021)] = silly_PID_xsec_constructor(0.76990E-02 * 1e3);
        result[iipair(2000004,1000021)] = silly_PID_xsec_constructor(0.36025E-02 * 1e3);
        result[iipair(2000005,1000021)] = silly_PID_xsec_constructor(0 * 1e3);
        //
        // Squark--anti-squark
        result[iipair(1000004,-1000004)] = silly_PID_xsec_constructor(0.395E-02 * 1e3);
        result[iipair(1000004,-1000003)] = silly_PID_xsec_constructor(0.150E-04 * 1e3);
        result[iipair(1000004,-1000001)] = silly_PID_xsec_constructor(0.139E-03 * 1e3);
        result[iipair(1000004,-1000002)] = silly_PID_xsec_constructor(0.380E-03 * 1e3);
        result[iipair(1000004,-2000002)] = silly_PID_xsec_constructor(0.171E-02 * 1e3);
        result[iipair(1000004,-2000001)] = silly_PID_xsec_constructor(0.681E-03 * 1e3);
        result[iipair(1000004,-2000003)] = silly_PID_xsec_constructor(0.804E-04 * 1e3);
        result[iipair(1000004,-2000004)] = silly_PID_xsec_constructor(0.375E-04 * 1e3);
        result[iipair(1000003,-1000003)] = silly_PID_xsec_constructor(0.391E-02 * 1e3);
        result[iipair(1000003,-1000001)] = silly_PID_xsec_constructor(0.234E-03 * 1e3);
        result[iipair(1000003,-1000002)] = silly_PID_xsec_constructor(0.611E-03 * 1e3);
        result[iipair(1000003,-2000002)] = silly_PID_xsec_constructor(0.273E-02 * 1e3);
        result[iipair(1000003,-2000001)] = silly_PID_xsec_constructor(0.113E-02 * 1e3);
        result[iipair(1000003,-2000003)] = silly_PID_xsec_constructor(0.156E-03 * 1e3);
        result[iipair(1000003,-2000004)] = silly_PID_xsec_constructor(0.786E-04 * 1e3);
        result[iipair(1000001,-1000001)] = silly_PID_xsec_constructor(0.497E-02 * 1e3);
        result[iipair(1000001,-1000002)] = silly_PID_xsec_constructor(0.159E-02 * 1e3);
        result[iipair(1000001,-2000002)] = silly_PID_xsec_constructor(0.697E-02 * 1e3);
        result[iipair(1000001,-2000001)] = silly_PID_xsec_constructor(0.404E-02 * 1e3);
        result[iipair(1000001,-2000003)] = silly_PID_xsec_constructor(0.113E-02 * 1e3);
        result[iipair(1000001,-2000004)] = silly_PID_xsec_constructor(0.668E-03 * 1e3);
        result[iipair(1000002,-1000002)] = silly_PID_xsec_constructor(0.660E-02 * 1e3);
        result[iipair(1000002,-2000002)] = silly_PID_xsec_constructor(0.877E-02 * 1e3);
        result[iipair(1000002,-2000001)] = silly_PID_xsec_constructor(0.710E-02 * 1e3);
        result[iipair(1000002,-2000003)] = silly_PID_xsec_constructor(0.277E-02 * 1e3);
        result[iipair(1000002,-2000004)] = silly_PID_xsec_constructor(0.171E-02 * 1e3);
        result[iipair(2000002,-2000002)] = silly_PID_xsec_constructor(0.835E-02 * 1e3);
        result[iipair(2000002,-2000001)] = silly_PID_xsec_constructor(0.196E-02 * 1e3);
        result[iipair(2000002,-2000003)] = silly_PID_xsec_constructor(0.753E-03 * 1e3);
        result[iipair(2000002,-2000004)] = silly_PID_xsec_constructor(0.462E-03 * 1e3);
        result[iipair(2000001,-2000001)] = silly_PID_xsec_constructor(0.665E-02 * 1e3);
        result[iipair(2000001,-2000003)] = silly_PID_xsec_constructor(0.297E-03 * 1e3);
        result[iipair(2000001,-2000004)] = silly_PID_xsec_constructor(0.174E-03 * 1e3);
        result[iipair(2000003,-2000003)] = silly_PID_xsec_constructor(0.531E-02 * 1e3);
        result[iipair(2000003,-2000004)] = silly_PID_xsec_constructor(0.192E-04 * 1e3);
        result[iipair(2000004,-2000004)] = silly_PID_xsec_constructor(0.515E-02 * 1e3);
        //
        result[iipair(1000005,-1000005)] = silly_PID_xsec_constructor(0.702E-02 * 1e3);
        result[iipair(1000006,-1000006)] = silly_PID_xsec_constructor(0.255E-01 * 1e3);
        result[iipair(2000005,-2000005)] = silly_PID_xsec_constructor(0.512E-02 * 1e3);
        result[iipair(2000006,-2000006)] = silly_PID_xsec_constructor(0.508E-02 * 1e3);
        //  
        // Squark--squark
        result[iipair(1000004,1000004)] = silly_PID_xsec_constructor(0.208E-04 * 1e3);
        result[iipair(1000004,1000003)] = silly_PID_xsec_constructor(0.122E-03 * 1e3);
        result[iipair(1000004,1000001)] = silly_PID_xsec_constructor(0.105E-02 * 1e3);
        result[iipair(1000004,1000002)] = silly_PID_xsec_constructor(0.270E-02 * 1e3);
        result[iipair(1000004,2000002)] = silly_PID_xsec_constructor(0.742E-03 * 1e3);
        result[iipair(1000004,2000001)] = silly_PID_xsec_constructor(0.278E-03 * 1e3);
        result[iipair(1000004,2000003)] = silly_PID_xsec_constructor(0.304E-04 * 1e3);
        result[iipair(1000004,2000004)] = silly_PID_xsec_constructor(0.132E-04 * 1e3);
        result[iipair(1000003,1000003)] = silly_PID_xsec_constructor(0.912E-04 * 1e3);
        result[iipair(1000003,1000001)] = silly_PID_xsec_constructor(0.215E-02 * 1e3);
        result[iipair(1000003,1000002)] = silly_PID_xsec_constructor(0.534E-02 * 1e3);
        result[iipair(1000003,2000002)] = silly_PID_xsec_constructor(0.160E-02 * 1e3);
        result[iipair(1000003,2000001)] = silly_PID_xsec_constructor(0.622E-03 * 1e3);
        result[iipair(1000003,2000003)] = silly_PID_xsec_constructor(0.665E-04 * 1e3);
        result[iipair(1000003,2000004)] = silly_PID_xsec_constructor(0.298E-04 * 1e3);
        result[iipair(1000001,1000001)] = silly_PID_xsec_constructor(0.759E-02 * 1e3);
        result[iipair(1000001,1000002)] = silly_PID_xsec_constructor(0.481E-01 * 1e3);
        result[iipair(1000001,2000002)] = silly_PID_xsec_constructor(0.165E-01 * 1e3);
        result[iipair(1000001,2000001)] = silly_PID_xsec_constructor(0.665E-02 * 1e3);
        result[iipair(1000001,2000003)] = silly_PID_xsec_constructor(0.622E-03 * 1e3);
        result[iipair(1000001,2000004)] = silly_PID_xsec_constructor(0.274E-03 * 1e3);
        result[iipair(1000002,1000002)] = silly_PID_xsec_constructor(0.410E-01 * 1e3);
        result[iipair(1000002,2000002)] = silly_PID_xsec_constructor(0.398E-01 * 1e3);
        result[iipair(1000002,2000001)] = silly_PID_xsec_constructor(0.166E-01 * 1e3);
        result[iipair(1000002,2000003)] = silly_PID_xsec_constructor(0.162E-02 * 1e3);
        result[iipair(1000002,2000004)] = silly_PID_xsec_constructor(0.742E-03 * 1e3);
        result[iipair(2000002,2000002)] = silly_PID_xsec_constructor(0.481E-01 * 1e3);
        result[iipair(2000002,2000001)] = silly_PID_xsec_constructor(0.582E-01 * 1e3);
        result[iipair(2000002,2000003)] = silly_PID_xsec_constructor(0.673E-02 * 1e3);
        result[iipair(2000002,2000004)] = silly_PID_xsec_constructor(0.341E-02 * 1e3);
        result[iipair(2000001,2000001)] = silly_PID_xsec_constructor(0.944E-02 * 1e3);
        result[iipair(2000001,2000003)] = silly_PID_xsec_constructor(0.279E-02 * 1e3);
        result[iipair(2000001,2000004)] = silly_PID_xsec_constructor(0.137E-02 * 1e3);
        result[iipair(2000003,2000003)] = silly_PID_xsec_constructor(0.122E-03 * 1e3);
        result[iipair(2000003,2000004)] = silly_PID_xsec_constructor(0.163E-03 * 1e3);
        result[iipair(2000004,2000004)] = silly_PID_xsec_constructor(0.276E-04 * 1e3);
        //
        // Associated production
        result[iipair(1000022,1000021)] = silly_PID_xsec_constructor(0.380E-02 * 1e3);
        result[iipair(1000023,1000021)] = silly_PID_xsec_constructor(0.444E-02 * 1e3);
        result[iipair(1000025,1000021)] = silly_PID_xsec_constructor(0.292E-04 * 1e3);
        result[iipair(1000035,1000021)] = silly_PID_xsec_constructor(0.111E-03 * 1e3);
        result[iipair(1000024,1000021)] = silly_PID_xsec_constructor(0.712E-02 * 1e3);
        result[iipair(-1000024,1000021)] = silly_PID_xsec_constructor(0.245E-02 * 1e3);
        result[iipair(1000037,1000021)] = silly_PID_xsec_constructor(0.214E-03 * 1e3);
        result[iipair(-1000037,1000021)] = silly_PID_xsec_constructor(0.706E-04 * 1e3);
        //
        result[iipair(1000022,1000004)] = silly_PID_xsec_constructor(0.132E-04 * 1e3);
        result[iipair(1000022,1000003)] = silly_PID_xsec_constructor(0.345E-04 * 1e3);
        result[iipair(1000022,1000001)] = silly_PID_xsec_constructor(0.206E-03 * 1e3);
        result[iipair(1000022,1000002)] = silly_PID_xsec_constructor(0.299E-03 * 1e3);
        result[iipair(1000022,2000002)] = silly_PID_xsec_constructor(0.707E-02 * 1e3);
        result[iipair(1000022,2000001)] = silly_PID_xsec_constructor(0.858E-03 * 1e3);
        result[iipair(1000022,2000003)] = silly_PID_xsec_constructor(0.147E-03 * 1e3);
        result[iipair(1000022,2000004)] = silly_PID_xsec_constructor(0.325E-03 * 1e3);
        result[iipair(1000023,1000004)] = silly_PID_xsec_constructor(0.301E-03 * 1e3);
        result[iipair(1000023,1000003)] = silly_PID_xsec_constructor(0.532E-03 * 1e3);
        result[iipair(1000023,1000001)] = silly_PID_xsec_constructor(0.340E-02 * 1e3);
        result[iipair(1000023,1000002)] = silly_PID_xsec_constructor(0.775E-02 * 1e3);
        result[iipair(1000023,2000002)] = silly_PID_xsec_constructor(0.688E-05 * 1e3);
        result[iipair(1000023,2000001)] = silly_PID_xsec_constructor(0.806E-06 * 1e3);
        result[iipair(1000023,2000003)] = silly_PID_xsec_constructor(0.129E-06 * 1e3);
        result[iipair(1000023,2000004)] = silly_PID_xsec_constructor(0.277E-06 * 1e3);
        result[iipair(1000025,1000004)] = silly_PID_xsec_constructor(0.258E-06 * 1e3);
        result[iipair(1000025,1000003)] = silly_PID_xsec_constructor(0.798E-06 * 1e3);
        result[iipair(1000025,1000001)] = silly_PID_xsec_constructor(0.565E-05 * 1e3);
        result[iipair(1000025,1000002)] = silly_PID_xsec_constructor(0.804E-05 * 1e3);
        result[iipair(1000025,2000002)] = silly_PID_xsec_constructor(0.293E-05 * 1e3);
        result[iipair(1000025,2000001)] = silly_PID_xsec_constructor(0.327E-06 * 1e3);
        result[iipair(1000025,2000003)] = silly_PID_xsec_constructor(0.470E-07 * 1e3);
        result[iipair(1000025,2000004)] = silly_PID_xsec_constructor(0.970E-07 * 1e3);
        result[iipair(1000035,1000004)] = silly_PID_xsec_constructor(0.662E-05 * 1e3);
        result[iipair(1000035,1000003)] = silly_PID_xsec_constructor(0.158E-04 * 1e3);
        result[iipair(1000035,1000001)] = silly_PID_xsec_constructor(0.113E-03 * 1e3);
        result[iipair(1000035,1000002)] = silly_PID_xsec_constructor(0.208E-03 * 1e3);
        result[iipair(1000035,2000002)] = silly_PID_xsec_constructor(0.148E-04 * 1e3);
        result[iipair(1000035,2000001)] = silly_PID_xsec_constructor(0.165E-05 * 1e3);
        result[iipair(1000035,2000003)] = silly_PID_xsec_constructor(0.235E-06 * 1e3);
        result[iipair(1000035,2000004)] = silly_PID_xsec_constructor(0.485E-06 * 1e3);
        result[iipair(1000024,1000004)] = silly_PID_xsec_constructor(0.475E-03 * 1e3);
        result[iipair(1000024,1000003)] = silly_PID_xsec_constructor(0.287E-03 * 1e3);
        result[iipair(1000024,1000001)] = silly_PID_xsec_constructor(0.140E-01 * 1e3);
        result[iipair(1000024,1000002)] = silly_PID_xsec_constructor(0.976E-03 * 1e3);
        result[iipair(1000024,2000002)] = silly_PID_xsec_constructor(0.00 * 1e3);
        result[iipair(1000024,2000001)] = silly_PID_xsec_constructor(0.00 * 1e3);
        result[iipair(1000024,2000003)] = silly_PID_xsec_constructor(0.00 * 1e3);
        result[iipair(1000024,2000004)] = silly_PID_xsec_constructor(0.00 * 1e3);
        result[iipair(1000037,1000004)] = silly_PID_xsec_constructor(0.749E-05 * 1e3);
        result[iipair(1000037,1000003)] = silly_PID_xsec_constructor(0.104E-04 * 1e3);
        result[iipair(1000037,1000001)] = silly_PID_xsec_constructor(0.629E-03 * 1e3);
        result[iipair(1000037,1000002)] = silly_PID_xsec_constructor(0.157E-04 * 1e3);
        result[iipair(1000037,2000002)] = silly_PID_xsec_constructor(0.00 * 1e3);
        result[iipair(1000037,2000001)] = silly_PID_xsec_constructor(0.00 * 1e3);
        result[iipair(1000037,2000003)] = silly_PID_xsec_constructor(0.00 * 1e3);
        result[iipair(1000037,2000004)] = silly_PID_xsec_constructor(0.00 * 1e3);
        result[iipair(-1000024,1000004)] = silly_PID_xsec_constructor(0.644E-03 * 1e3);
        result[iipair(-1000024,1000003)] = silly_PID_xsec_constructor(0.287E-03 * 1e3);
        result[iipair(-1000024,1000001)] = silly_PID_xsec_constructor(0.775E-03 * 1e3);
        result[iipair(-1000024,1000002)] = silly_PID_xsec_constructor(0.618E-02 * 1e3);
        result[iipair(-1000024,2000002)] = silly_PID_xsec_constructor(0.00 * 1e3);
        result[iipair(-1000024,2000001)] = silly_PID_xsec_constructor(0.00 * 1e3);
        result[iipair(-1000024,2000003)] = silly_PID_xsec_constructor(0.00 * 1e3);
        result[iipair(-1000024,2000004)] = silly_PID_xsec_constructor(0.00 * 1e3);
        result[iipair(-1000037,1000004)] = silly_PID_xsec_constructor(0.109E-04 * 1e3);
        result[iipair(-1000037,1000003)] = silly_PID_xsec_constructor(0.104E-04 * 1e3);
        result[iipair(-1000037,1000001)] = silly_PID_xsec_constructor(0.294E-04 * 1e3);
        result[iipair(-1000037,1000002)] = silly_PID_xsec_constructor(0.115E-03 * 1e3);
        result[iipair(-1000037,2000002)] = silly_PID_xsec_constructor(0.00 * 1e3);
        result[iipair(-1000037,2000001)] = silly_PID_xsec_constructor(0.00 * 1e3);
        result[iipair(-1000037,2000003)] = silly_PID_xsec_constructor(0.00 * 1e3);
        result[iipair(-1000037,2000004)] = silly_PID_xsec_constructor(0.00 * 1e3);
        //
        // EW
        result[iipair(1000022,1000022)] = silly_PID_xsec_constructor(0.990E-03 * 1e3);
        result[iipair(1000022,1000023)] = silly_PID_xsec_constructor(0.137E-03 * 1e3);
        result[iipair(1000022,1000025)] = silly_PID_xsec_constructor(0.744E-04 * 1e3);
        result[iipair(1000022,1000035)] = silly_PID_xsec_constructor(0.167E-04 * 1e3);
        result[iipair(1000022,1000024)] = silly_PID_xsec_constructor(0.350E-03 * 1e3);
        result[iipair(1000022,1000037)] = silly_PID_xsec_constructor(0.107E-03 * 1e3);
        result[iipair(1000022,-1000024)] = silly_PID_xsec_constructor(0.162E-03 * 1e3);
        result[iipair(1000022,-1000037)] = silly_PID_xsec_constructor(0.460E-04 * 1e3);
        result[iipair(1000023,1000023)] = silly_PID_xsec_constructor(0.142E-02 * 1e3);
        result[iipair(1000023,1000025)] = silly_PID_xsec_constructor(0.282E-03 * 1e3);
        result[iipair(1000023,1000035)] = silly_PID_xsec_constructor(0.678E-04 * 1e3);
        result[iipair(1000023,1000024)] = silly_PID_xsec_constructor(0.411E-01 * 1e3);
        result[iipair(1000023,1000037)] = silly_PID_xsec_constructor(0.183E-04 * 1e3);
        result[iipair(1000023,-1000024)] = silly_PID_xsec_constructor(0.191E-01 * 1e3);
        result[iipair(1000023,-1000037)] = silly_PID_xsec_constructor(0.719E-05 * 1e3);
        result[iipair(1000025,000025)] = silly_PID_xsec_constructor(0.141E-07 * 1e3);
        result[iipair(1000025,1000035)] = silly_PID_xsec_constructor(0.175E-02 * 1e3);
        result[iipair(1000025,1000024)] = silly_PID_xsec_constructor(0.412E-03 * 1e3);
        result[iipair(1000025,1000037)] = silly_PID_xsec_constructor(0.256E-02 * 1e3);
        result[iipair(1000025,-1000024)] = silly_PID_xsec_constructor(0.169E-03 * 1e3);
        result[iipair(1000025,-1000037)] = silly_PID_xsec_constructor(0.968E-03 * 1e3);
        result[iipair(1000035,1000035)] = silly_PID_xsec_constructor(0.834E-06 * 1e3);
        result[iipair(1000035,1000024)] = silly_PID_xsec_constructor(0.248E-04 * 1e3);
        result[iipair(1000035,1000037)] = silly_PID_xsec_constructor(0.248E-02 * 1e3);
        result[iipair(1000035,-1000024)] = silly_PID_xsec_constructor(0.968E-05 * 1e3);
        result[iipair(1000035,-1000037)] = silly_PID_xsec_constructor(0.935E-03 * 1e3);
        result[iipair(1000024,-1000024)] = silly_PID_xsec_constructor(0.322E-01 * 1e3);
        result[iipair(1000024,-1000037)] = silly_PID_xsec_constructor(0.130E-03 * 1e3);
        result[iipair(1000037,-1000024)] = silly_PID_xsec_constructor(0.130E-03 * 1e3);
        result[iipair(1000037,-1000037)] = silly_PID_xsec_constructor(0.183E-02 * 1e3);
        //
        // Sleptons TODO: What do to with first and second generation equals
        result[iipair(1000011,-1000011)] = silly_PID_xsec_constructor(0.210E-02 * 1e3); // Equal for ~e and ~mu in Prospino
        result[iipair(1000013,-1000013)] = silly_PID_xsec_constructor(0.210E-02 * 1e3);
        result[iipair(2000011,-2000011)] = silly_PID_xsec_constructor(0.502E-02 * 1e3); // Equal for ~e and ~mu in Prospino
        result[iipair(2000013,-2000013)] = silly_PID_xsec_constructor(0.502E-02 * 1e3);
        result[iipair(1000012,-1000012)] = silly_PID_xsec_constructor(0.217E-02 * 1e3); // Equal for ~nu_e and ~nu_mu in Prospino
        result[iipair(1000014,-1000014)] = silly_PID_xsec_constructor(0.217E-02 * 1e3);
        result[iipair(-1000011,1000012)] = silly_PID_xsec_constructor(0.553E-02 * 1e3); // Equal for ~nu_e and ~nu_mu in Prospino
        result[iipair(-1000013,1000014)] = silly_PID_xsec_constructor(0.553E-02 * 1e3);
        result[iipair(1000011,-1000012)] = silly_PID_xsec_constructor(0.245E-02 * 1e3);
        result[iipair(1000013,-1000014)] = silly_PID_xsec_constructor(0.245E-02 * 1e3);
        result[iipair(1000015,-1000015)] = silly_PID_xsec_constructor(0.556E-02 * 1e3);
        result[iipair(2000015,-2000015)] = silly_PID_xsec_constructor(0.201E-02 * 1e3);
        result[iipair(1000015,-2000015)] = silly_PID_xsec_constructor(0.907E-04 * 1e3); // Equal for ~tau_1^- ~tau_2^+ and its cc in Prospino. Suspicious difference in value compared to Pythia LO
        result[iipair(-1000015,2000015)] = silly_PID_xsec_constructor(0.907E-04 * 1e3);
        result[iipair(1000016,-1000016)] = silly_PID_xsec_constructor(0.220E-02 * 1e3);
        result[iipair(-1000015,1000016)] = silly_PID_xsec_constructor(0.262E-03 * 1e3);
        result[iipair(1000015,-1000016)] = silly_PID_xsec_constructor(0.124E-03 * 1e3);
        result[iipair(-2000015,1000016)] = silly_PID_xsec_constructor(0.542E-02 * 1e3);
      }

    }



    /// Get a PIDPairCrossSectionFunc pointing to PIDPairCrossSection_dummy,
    /// with the second argument already filled by *Dep::MSSM_spectrum
    void getPIDPairCrossSectionFunc_dummy(PIDPairCrossSectionFuncType& result)
    {
      using namespace Pipes::getPIDPairCrossSectionFunc_dummy;
      result = std::bind(PIDPairCrossSection_dummy, std::placeholders::_1, *Dep::MSSM_spectrum);
    }





    /// Get the cross-section from the xsec_example backend
    xsec_container PIDPairCrossSection_xsec_example(iipair pids, 
                                          const Spectrum& MSSM_spectrum)
                                          // double (*xsec_fb_fptr)(iipair&, map_str_dbl&, map_str_bool&))
    {
      xsec_container xs_result;

      // Get an SLHA1 object
      const SLHAstruct& slha = MSSM_spectrum.getSLHAea(1);

      // Call cross-section calculation 
      map_str_dbl proc_params;
      map_str_bool proc_flags;

      double xs_fb = 3.1415;

      // double xs_fb = 1.0e-3 * xsec_fb_fptr(pids, proc_params, proc_flags);

      double xs_rel_err = 0.01;

      // Save result in xs_result
      xs_result.set_xsec(xs_fb, xs_rel_err);

      // Construct info string of the form "PID1:<PID1>, PID2:<PID2>"
      std::stringstream info_ss;
      info_ss << "PID1:" << pids.first << ", " << "PID2:" << pids.second;
      xs_result.set_info_string(info_ss.str());

      return xs_result;
    }


    /// Get a PIDPairCrossSectionFunc pointing to the xsec_example backend
    void getPIDPairCrossSectionFunc_xsec_example(PIDPairCrossSectionFuncType& result)
    {
      using namespace Pipes::getPIDPairCrossSectionFunc_xsec_example;

      // result = std::bind(PIDPairCrossSection_xsec_example, std::placeholders::_1, *Dep::MSSM_spectrum, BEreq::xsec_example_xsec_fb.pointer());
      result = std::bind(PIDPairCrossSection_xsec_example, std::placeholders::_1, *Dep::MSSM_spectrum);

      // Bind to the backend function xsec_fb with the two first arguments fixed

      // // Get an SLHA1 object
      // const SLHAstruct& slha = MSSM_spectrum.getSLHAea(1);

      // _Anders
      // BACKEND_REQ(xsec_example_xsec_fb, (), double, (std::vector<double>&, map_str_dbl&, map_str_bool&))

      // map_str_dbl proc_params;
      // map_str_bool proc_flags;

      // result = std::bind(BEreq::xsec_example_xsec_fb.pointer(), std::placeholders::_1, proc_params, proc_flags);
    }





    /// Get a map between Pythia process codes and cross-sections
    void getProcessCrossSectionsMap(map_int_process_xsec& result)
    {
      using namespace Pipes::getProcessCrossSectionsMap;

      // Use a static variable to communicate the result calculated on thread 0 during 
      // iteration XSEC_CALCULATION to all threads during iteration START_SUBPROCESS
      static map_int_process_xsec shared_result;

      const static bool set_missing_xsecs_to_zero = runOptions->getValueOrDef<bool>(false, "set_missing_xsecs_to_zero");

      // Only thread 0
      if(*Loop::iteration == COLLIDER_INIT)
      {
        shared_result.clear();
      }

      // All threads
      if(*Loop::iteration == COLLIDER_INIT_OMP)
      {
        result.clear();
      }

      // Only thread 0
      if(*Loop::iteration == XSEC_CALCULATION)
      {
        cout << DEBUG_PREFIX << "getProcessCrossSectionsMap: it = XSEC_CALCULATION, ActiveProcessCodes.size() = " << Dep::ActiveProcessCodes->size() << endl;          
        cout << DEBUG_PREFIX << "getProcessCrossSectionsMap: it = XSEC_CALCULATION, ActiveProcessCodeToPIDPairsMap.size() = " << Dep::ActiveProcessCodeToPIDPairsMap->size() << endl;          

        // Loop over all active processes and construct the cross-section map (shared_result)
        for (size_t i = 0; i != Dep::ActiveProcessCodes->size(); ++i)
        {
          // Get process code
          int proc_code = Dep::ActiveProcessCodes->at(i);

          // Construct a process_xsec_container instance to be stored in the shared_result map
          process_xsec_container proc_xs;
          proc_xs.set_process_code(proc_code);

          // Get iterator bounds (as a pair) over the multimap entries that match the key proc_code
          auto mm_range = Dep::ActiveProcessCodeToPIDPairsMap->equal_range(proc_code);

          // Loop over these elements in the multimap
          for (auto mm_it = mm_range.first; mm_it != mm_range.second; ++mm_it)
          {
            const iipair& pids = mm_it->second;

            // Obtain the cross-section from the PID pair via the PIDPairCrossSectionsMap (map_iipair_PID_pair_xsec) dependency
            cout << DEBUG_PREFIX << "Looking up PID pair: " << pids.first << "," << pids.second << endl;

            PID_pair_xsec_container pids_xs;
            map_iipair_PID_pair_xsec::const_iterator iter = Dep::PIDPairCrossSectionsMap->find(pids);
            if (iter != Dep::PIDPairCrossSectionsMap->end())
            {
              pids_xs = iter->second;
              cout << DEBUG_PREFIX << "--> Got it! Cross-section is: " << pids_xs.xsec() << " fb" << endl;
            }
            else
            {
              if(set_missing_xsecs_to_zero)
              {
                cout << DEBUG_PREFIX << "--> Not found! Creating 0-valued xsec_container " << endl;
                pids_xs.set_xsec(0.0, 0.0);
              }
              else
              {
                std::stringstream errmsg_ss;
                errmsg_ss << "No cross-section provided for PID pair [" << pids.first << "," << pids.second <<"]. ";
                ColliderBit_error().raise(LOCAL_INFO, errmsg_ss.str());
              }
            }

            // Accumulate result in the process_xsec_container proc_xs
            proc_xs.sum_xsecs(pids_xs.xsec(), pids_xs.xsec_err());
            proc_xs.add_related_PID_pair(pids);
          }

          // Construct info string of the form "ProcessCode:<proc_code>"
          std::stringstream info_ss;
          info_ss << "ProcessCode:" << proc_code;
          proc_xs.set_info_string(info_ss.str());

          // Now we figure out if the proc_code process shares the cross-section
          // stored in in proc_xs with any other process codes

          // Loop over *all* entries (process code <--> PID pair) in the multimap all_process_codes_to_PID_pairs
          for (auto mm_it = all_process_codes_to_PID_pairs.begin(); mm_it != all_process_codes_to_PID_pairs.end(); ++mm_it)
          {
            // Extract the process code (pc) and PID pair (pp)
            int pc = mm_it->first;
            const iipair& pp = mm_it->second;

            if (pc == proc_code) continue;

            // @todo What's the right choice here?
            // // Don't add more copies of the same process code! ...Or should we?
            // if(std::find(proc_xs.processes_sharing_xsec().begin(), proc_xs.processes_sharing_xsec().end(), pc) != proc_xs.processes_sharing_xsec().end()) 

            // Check if the PID pair pp mathces one of the PID pairs for the proc_code process
            if(std::find(proc_xs.related_PID_pairs().begin(), proc_xs.related_PID_pairs().end(), pp) != proc_xs.related_PID_pairs().end()) 
            {
              // Check that pc is itself in one of the active processes, i.e. listed in Dep::ActiveProcessCodes
              if(std::find(Dep::ActiveProcessCodes->begin(), Dep::ActiveProcessCodes->end(), pc) != Dep::ActiveProcessCodes->end())  
              {
                // Add pc to the list of processes that share cross-section with proc_code
                proc_xs.add_process_sharing_xsec(pc);
              }
              else
              {
                std::stringstream errmsg_ss;
                errmsg_ss << "For correct cross-section scaling of collider process " << proc_code;
                errmsg_ss << ", process " << pc << " must also be activated. Please check your collider settings.";
                ColliderBit_error().raise(LOCAL_INFO, errmsg_ss.str());
              }
            }
          }

          // Store proc_xs in the shared_result map
          shared_result[proc_code] = proc_xs;
        }

        // Let thread 0 return the correct result already after iteration XSEC_CALCULATION
        result = shared_result;
      }


      // All threads
      if (*Loop::iteration == START_SUBPROCESS)
      {
        // All threads read the result from shared_result
        result = shared_result;
      }

    }


    /// Compute a cross-section from Monte Carlo
    void getMCxsec(MC_xsec_container& result)
    {
      using namespace Pipes::getMCxsec;

      // Don't bother if there are no analyses that will use this.
      if (Dep::RunMC->analyses.empty()) return;

      // Reset the xsec objects on all threads
      if (*Loop::iteration == COLLIDER_INIT_OMP)
      {
        result.reset();
      }

      // If we are in the main event loop, count the event towards cross-section normalisation on this thread
      if (*Loop::iteration >= 0)
      {
        result.log_event();
      }

      // Extract the xsecs from the MC on each thread
      if (*Loop::iteration == END_SUBPROCESS && Dep::RunMC->event_generation_began)
      {
        if (not Dep::RunMC->exceeded_maxFailedEvents)
        {
          const double xs_fb = (*Dep::HardScatteringSim)->xsec_fb();
          const double xserr_fb = (*Dep::HardScatteringSim)->xsecErr_fb();
          result.set_xsec(xs_fb, xserr_fb);
          #ifdef COLLIDERBIT_DEBUG
            cout << DEBUG_PREFIX << "xs_fb = " << xs_fb << " +/- " << xserr_fb << endl;
          #endif
        }
      }

      // Gather the xsecs from all threads into one
      if (*Loop::iteration == COLLIDER_FINALIZE)
      {
        result.gather_xsecs();
      }

    }

    /// Return MC_xsec_container as const base_xsec_container*
    void get_MC_xsec_as_base(const base_xsec_container*& result)
    {
      using namespace Pipes::get_MC_xsec_as_base;
      result = &(*Dep::TotalCrossSectionFromMC);
    }

    /// Return xsec_container as const base_xsec_container*
    void get_xsec_as_base(const base_xsec_container*& result)
    {
      using namespace Pipes::get_xsec_as_base;
      result = &(*Dep::TotalCrossSection);
    }



    /// Get a cross-section from NLL-FAST
    void getNLLFastxsec(xsec_container& result)
    {
      using namespace Pipes::getNLLFastxsec;

      // Use a static variable to communicate the result calculated on thread 0 during 
      // iteration XSEC_CALCULATION to all threads during iteration START_SUBPROCESS
      static xsec_container shared_result;

      // Don't bother if there are no analyses that will use this.
      if (Dep::RunMC->analyses.empty()) return;

      // Only thread 0
      if(*Loop::iteration == COLLIDER_INIT) shared_result.reset();
      
      // All threads
      if (*Loop::iteration == COLLIDER_INIT_OMP) result.reset();

      // Set the xsec and its error.
      // Only thread 0
      if (*Loop::iteration == XSEC_CALCULATION)
      {
        double xs_fb = 0.1;             // replace with xsec from NLL-Fast
        double xserr_fb = 0.1 * xs_fb;  // or whatever
        shared_result.set_xsec(xs_fb, xserr_fb);

        // Let thread 0 return the correct result already after iteration XSEC_CALCULATION
        result = shared_result;
      }

      // All threads
      if (*Loop::iteration == START_SUBPROCESS)
      {
        // All threads copy the result from shared_result
        result = shared_result;
      }

    }


    /// A function that reads the total cross-section from the input file, but builds up the number of events from the event loop
    void getYAMLxsec(xsec_container& result)
    {
      using namespace Pipes::getYAMLxsec;

      // Use a static variable to communicate the result calculated on thread 0 during 
      // iteration XSEC_CALCULATION to all threads during iteration START_SUBPROCESS
      static xsec_container shared_result;

      // Don't bother if there are no analyses that will use this.
      if (Dep::RunMC->analyses.empty()) return;

      static std::pair<str,str> xsec_pnames;
      static str input_unit; 
      static bool input_fractional_uncert = false;

      static bool first = true;
      if (*Loop::iteration == BASE_INIT)
      {

        if (first)
        {
          // Determine the correct combination of parameters
          if ((runOptions->hasKey("xsec_fb")) && (runOptions->hasKey("xsec_uncert_fb")))
          {
            xsec_pnames.first = "xsec_fb";
            xsec_pnames.second = "xsec_uncert_fb";
            input_unit = "fb";
            input_fractional_uncert = false;
          }
          else if ((runOptions->hasKey("xsec_fb")) && (runOptions->hasKey("xsec_fractional_uncert")))
          {
            xsec_pnames.first = "xsec_fb";
            xsec_pnames.second = "xsec_fractional_uncert";
            input_unit = "fb";
            input_fractional_uncert = true;
          }
          else if ((runOptions->hasKey("xsec_pb")) && (runOptions->hasKey("xsec_uncert_pb")))
          {
            xsec_pnames.first = "xsec_pb";
            xsec_pnames.second = "xsec_uncert_pb";
            input_unit = "pb";
            input_fractional_uncert = false;
          }
          else if ((runOptions->hasKey("xsec_pb")) && (runOptions->hasKey("xsec_fractional_uncert")))
          {
            xsec_pnames.first = "xsec_pb";
            xsec_pnames.second = "xsec_fractional_uncert";
            input_unit = "pb";
            input_fractional_uncert = true;
          }
          else
          {
            std::stringstream errmsg_ss;
            errmsg_ss << "Unknown combination of options for function getYAMLxsec." << endl;
            errmsg_ss << "Needs one of the following sets of option names:" << endl;
            errmsg_ss << "  xsec_fb, xsec_uncert_fb" << endl;
            errmsg_ss << "  xsec_fb, xsec_fractional_uncert" << endl;
            errmsg_ss << "  xsec_pb, xsec_uncert_pb" << endl;
            errmsg_ss << "  xsec_pb, xsec_fractional_uncert" << endl;
            ColliderBit_error().raise(LOCAL_INFO, errmsg_ss.str());
          }

          first = false;
        }
      }

      // Retrieve the total cross-section and cross-section error
      const static double input_xsec = runOptions->getValue<double>(xsec_pnames.first);
      const static double input_xsec_uncert = runOptions->getValue<double>(xsec_pnames.second);

      // Only thread 0
      if(*Loop::iteration == COLLIDER_INIT) shared_result.reset();
      
      // All threads
      if (*Loop::iteration == COLLIDER_INIT_OMP) result.reset();

      // Set the xsec and its error
      // Only thread 0
      if (*Loop::iteration == XSEC_CALCULATION)
      {
        std::pair<double,double> temp = convert_xsecs_to_fb(input_xsec, input_xsec_uncert, input_unit, input_fractional_uncert);
        double xsec_fb = temp.first;
        double xsec_uncert_fb = temp.second;
        shared_result.set_xsec(xsec_fb, xsec_uncert_fb);

        // Let thread 0 return the correct result already after iteration XSEC_CALCULATION
        result = shared_result;
      }

      // All threads
      if (*Loop::iteration == START_SUBPROCESS)
      {
        // All threads copy the result from shared_result
        result = shared_result;
      }

    }



    /// A function that reads a list of (SLHA file, total cross-section) pairs from the input YAML file
    void getYAMLxsec_SLHA(xsec_container& result)
    {
      using namespace Pipes::getYAMLxsec_SLHA;

      // Use a static variable to communicate the result calculated on thread 0 during 
      // iteration XSEC_CALCULATION to all threads during iteration START_SUBPROCESS
      static xsec_container shared_result;

      // Don't bother if there are no analyses that will use this.
      if (Dep::RunMC->analyses.empty()) return;

      static std::pair<str,str> xsec_pnames;
      static str input_unit; 
      static bool input_fractional_uncert = false;

      static bool first = true;
      if (*Loop::iteration == BASE_INIT)
      {

        if (first)
        {
          // Determine the correct combination of parameters
          if ((runOptions->hasKey("xsec_fb")) && (runOptions->hasKey("xsec_uncert_fb")))
          {
            xsec_pnames.first = "xsec_fb";
            xsec_pnames.second = "xsec_uncert_fb";
            input_unit = "fb";
            input_fractional_uncert = false;
          }
          else if ((runOptions->hasKey("xsec_fb")) && (runOptions->hasKey("xsec_fractional_uncert")))
          {
            xsec_pnames.first = "xsec_fb";
            xsec_pnames.second = "xsec_fractional_uncert";
            input_unit = "fb";
            input_fractional_uncert = true;
          }
          else if ((runOptions->hasKey("xsec_pb")) && (runOptions->hasKey("xsec_uncert_pb")))
          {
            xsec_pnames.first = "xsec_pb";
            xsec_pnames.second = "xsec_uncert_pb";
            input_unit = "pb";
            input_fractional_uncert = false;
          }
          else if ((runOptions->hasKey("xsec_pb")) && (runOptions->hasKey("xsec_fractional_uncert")))
          {
            xsec_pnames.first = "xsec_pb";
            xsec_pnames.second = "xsec_fractional_uncert";
            input_unit = "pb";
            input_fractional_uncert = true;
          }
          else
          {
            std::stringstream errmsg_ss;
            errmsg_ss << "Unknown combination of options for function getYAMLxsec_SLHA." << endl;
            errmsg_ss << "Needs one of the following sets of option names:" << endl;
            errmsg_ss << "  xsec_fb, xsec_uncert_fb" << endl;
            errmsg_ss << "  xsec_fb, xsec_fractional_uncert" << endl;
            errmsg_ss << "  xsec_pb, xsec_uncert_pb" << endl;
            errmsg_ss << "  xsec_pb, xsec_fractional_uncert" << endl;
            ColliderBit_error().raise(LOCAL_INFO, errmsg_ss.str());
          }

          first = false;
        }
      }

      // Get the filename, look for it in the xsec and uncertainty lists
      const static YAML::Node colNode_xsec = runOptions->getValue<YAML::Node>(xsec_pnames.first);
      const static Options colOptions_xsec(colNode_xsec);
      const static YAML::Node colNode_uncert = runOptions->getValue<YAML::Node>(xsec_pnames.second);
      const static Options colOptions_uncert(colNode_uncert);
      static str filename;

      if (*Loop::iteration == BASE_INIT)
      {
        // Update the SLHA filename
        filename = Dep::SLHAFileNameAndContent->first;

        // Look for the filename in the xsec lists
        if (!colOptions_xsec.hasKey(filename)) piped_invalid_point.request(str("No cross-section found for SLHA file ").append(filename));
        if (!colOptions_uncert.hasKey(filename)) piped_invalid_point.request(str("No fractional cross-section uncertainty found for SLHA file ").append(filename));
      }

      // Only thread 0
      if(*Loop::iteration == COLLIDER_INIT) shared_result.reset();
      
      // All threads
      if (*Loop::iteration == COLLIDER_INIT_OMP) result.reset();

      // Set the xsec and its error
      if (*Loop::iteration == XSEC_CALCULATION)
      {
        double input_xsec = colOptions_xsec.getValue<double>(filename);
        double input_xsec_uncert = colOptions_uncert.getValue<double>(filename);

        std::pair<double,double> temp = convert_xsecs_to_fb(input_xsec, input_xsec_uncert, input_unit, input_fractional_uncert);
        double xsec_fb = temp.first;
        double xsec_uncert_fb = temp.second;
        shared_result.set_xsec(xsec_fb, xsec_uncert_fb);

        // Let thread 0 return the correct result already after iteration XSEC_CALCULATION
        result = shared_result;
      }

      // All threads
      if (*Loop::iteration == START_SUBPROCESS)
      {
        // All threads copy the result from shared_result
        result = shared_result;
      }

    }


    /// A function that assigns a total cross-sections directly from the scan parameters
    /// (for models CB_SLHA_simpmod_scan_model and CB_SLHA_scan_model)
    void getYAMLxsec_param(xsec_container& result)
    {
      using namespace Pipes::getYAMLxsec_param;

      // Use a static variable to communicate the result calculated on thread 0 during 
      // iteration XSEC_CALCULATION to all threads during iteration START_SUBPROCESS
      static xsec_container shared_result;

      // Don't bother if there are no analyses that will use this.
      if (Dep::RunMC->analyses.empty()) return;

      static std::vector<str> pnames;
      static std::pair<str,str> xsec_pnames;

      static str input_unit; 
      static bool input_fractional_uncert = false;

      static bool first = true;
      if (*Loop::iteration == BASE_INIT)
      {

        if (first)
        {

          // Get all parameter names
          for (const auto& parname_parptr_pair : Param)
          {
            pnames.push_back(parname_parptr_pair.first);
          }

          // Determine the correct combination of parameters
          if ((std::find(pnames.begin(), pnames.end(), "xsec_fb") != pnames.end()) 
               && (std::find(pnames.begin(), pnames.end(), "xsec_uncert_fb") != pnames.end()))
          {
            xsec_pnames.first = "xsec_fb";
            xsec_pnames.second = "xsec_uncert_fb";
            input_unit = "fb";
            input_fractional_uncert = false;
          }
          else if ((std::find(pnames.begin(), pnames.end(), "xsec_fb") != pnames.end()) 
                    && (std::find(pnames.begin(), pnames.end(), "xsec_fractional_uncert") != pnames.end()))
          {
            xsec_pnames.first = "xsec_fb";
            xsec_pnames.second = "xsec_fractional_uncert";
            input_unit = "fb";
            input_fractional_uncert = true;
          }
          else if ((std::find(pnames.begin(), pnames.end(), "xsec_pb") != pnames.end()) 
                    && (std::find(pnames.begin(), pnames.end(), "xsec_uncert_pb") != pnames.end()))
          {
            xsec_pnames.first = "xsec_pb";
            xsec_pnames.second = "xsec_uncert_pb";
            input_unit = "pb";
            input_fractional_uncert = false;
          }
          else if ((std::find(pnames.begin(), pnames.end(), "xsec_pb") != pnames.end()) 
                    && (std::find(pnames.begin(), pnames.end(), "xsec_fractional_uncert") != pnames.end()))
          {
            xsec_pnames.first = "xsec_pb";
            xsec_pnames.second = "xsec_fractional_uncert";
            input_unit = "pb";
            input_fractional_uncert = true;
          }
          else
          {
            std::stringstream errmsg_ss;
            errmsg_ss << "Unknown combination of parameters for function getYAMLxsec_param." << endl;
            errmsg_ss << "Needs one of the following sets of parameter names:" << endl;
            errmsg_ss << "  xsec_fb, xsec_uncert_fb" << endl;
            errmsg_ss << "  xsec_fb, xsec_fractional_uncert" << endl;
            errmsg_ss << "  xsec_pb, xsec_uncert_pb" << endl;
            errmsg_ss << "  xsec_pb, xsec_fractional_uncert" << endl;
            ColliderBit_error().raise(LOCAL_INFO, errmsg_ss.str());
          }

          first = false;
        }
      }

      // Only thread 0
      if(*Loop::iteration == COLLIDER_INIT) shared_result.reset();
      
      // All threads
      if (*Loop::iteration == COLLIDER_INIT_OMP) result.reset();

      // Set the xsec and its error
      // Only thread 0
      if (*Loop::iteration == XSEC_CALCULATION)
      {
        double input_xsec = *Param.at(xsec_pnames.first);
        double input_xsec_uncert = *Param.at(xsec_pnames.second); 

        std::pair<double,double> temp = convert_xsecs_to_fb(input_xsec, input_xsec_uncert, input_unit, input_fractional_uncert);
        double xsec_fb = temp.first;
        double xsec_uncert_fb = temp.second;
        shared_result.set_xsec(xsec_fb, xsec_uncert_fb);

        // Let thread 0 return the correct result already after iteration XSEC_CALCULATION
        result = shared_result;
      }

      // All threads
      if (*Loop::iteration == START_SUBPROCESS)
      {
        // All threads copy the result from shared_result
        result = shared_result;
      }

    }


    /// Get cross-section info as map_str_dbl (for simple printing)
    void getTotalCrossSectionAsMap(map_str_dbl& result)
    {
      using namespace Pipes::getTotalCrossSectionAsMap;

      // @todo Do we need this to ensure that the result map is always of the same length (for the printer)?
      // // Append the xsec info for the current collider to the result map
      // if (*Loop::iteration == COLLIDER_INIT)
      // {
      //   xsec empty_xs;
      //   for(auto s_d_pair : empty_xs.get_content_as_map())
      //   {
      //     std::string new_key(Dep::RunMC->current_collider());
      //     new_key.append("__").append(s_d_pair.first);
      //     result[new_key] = s_d_pair.second;
      //   }
      // }

      // Append the xsec info for the current collider to the result map
      if (*Loop::iteration == COLLIDER_FINALIZE)
      {
        for(auto s_d_pair : (*Dep::TotalCrossSection)->get_content_as_map())
        {
          std::string new_key(Dep::RunMC->current_collider());
          new_key.append("__").append(s_d_pair.first);
          result[new_key] = s_d_pair.second;
        }
      }
    }

  }
}
