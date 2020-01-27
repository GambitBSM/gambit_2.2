//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  Functions of module FlavBit
///
///  *********************************************
///
///  Authors (add name and date if you modify):
///
///  \author Nazila Mahmoudi
///  \date 2013 Oct
///  \date 2014
///  \date 2015 Feb
///  \date 2016 Jul
///  \date 2018 Jan
///  \date 2019 Aug
///
///  \author Marcin Chrzaszcz
///  \date 2015 May
///  \date 2015 July
///  \date 2015 August
///  \date 2016 July
///  \date 2016 August
///  \date 2016 October
///  \date 2018 Jan
///
///  \author Anders Kvellestad
///          (anders.kvellestad@fys.uio.no)
///  \date 2013 Nov
///
///  \author Pat Scott
///          (pat.scott@uq.edu.au)
///  \date 2015 May, June
///  \date 2016 Aug
///  \date 2017 March
///  \date 2019 Oct
///
///  \author Tomas Gonzalo
///          (t.e.gonzalo@fys.uio.no)
///  \date 2017 July
///
///  \date 2017 July
///  \author Jihyun Bhom
///          (jihyun.bhom@ifj.edu.pl)
///  \date 2019 July
///  \date 2019 Nov
///  \date 2019 Dec
///
///  \author Markus Prim
///          (markus.prim@kit.edu)
///  \date 2019 August
///
///  *********************************************

#include <string>
#include <iostream>
#include <fstream>
#include <map>

#include "gambit/Elements/gambit_module_headers.hpp"
#include "gambit/FlavBit/FlavBit_rollcall.hpp"
#include "gambit/FlavBit/FlavBit_types.hpp"
#include "gambit/FlavBit/Flav_reader.hpp"
#include "gambit/FlavBit/Kstarmumu_theory_err.hpp"
#include "gambit/FlavBit/flav_utils.hpp"
#include "gambit/FlavBit/flav_loop_functions.hpp"
#include "gambit/Elements/spectrum.hpp"
#include "gambit/Utils/statistics.hpp"
#include "gambit/cmake/cmake_variables.hpp"


#define FLAVBIT_DEBUG
//#define FLAVBIT_DEBUG_LL

namespace Gambit
{

  namespace FlavBit
  {

    using namespace std;
    namespace ublas = boost::numeric::ublas;

    const bool flav_debug =
    #ifdef FLAVBIT_DEBUG
      true;
    #else
      false;
    #endif

    const bool flav_debug_LL =
    #ifdef FLAVBIT_DEBUG_LL
      true;
    #else
      false;
    #endif

    /// Find the path to the latest installed version of the HepLike data
    str path_to_latest_heplike_data()
    {
      std::vector<str> working_data = Backends::backendInfo().working_versions("HepLikeData");
      if (working_data.empty()) FlavBit_error().raise(LOCAL_INFO, "No working HepLikeData installations detected.");
      std::sort(working_data.begin(), working_data.end());
      return Backends::backendInfo().corrected_path("HepLikeData", working_data.back());
    }

    /// Fill SuperIso model info structure
    void SI_fill(parameters &result)
    {
      using namespace Pipes::SI_fill;
      using namespace std;

      SLHAstruct spectrum;
      // Obtain SLHAea object from spectrum
      if (ModelInUse("WC"))
      {
        spectrum = Dep::SM_spectrum->getSLHAea(2);
      }
      else if (ModelInUse("MSSM63atMGUT") or ModelInUse("MSSM63atQ"))
      {
        spectrum = Dep::MSSM_spectrum->getSLHAea(2);
        // Add the MODSEL block if it is not provided by the spectrum object.
        SLHAea_add(spectrum,"MODSEL",1, 0, "General MSSM", false);
      }
      else
      {
        FlavBit_error().raise(LOCAL_INFO, "Unrecognised model.");
      }

      BEreq::Init_param(&result);

      int ie,je;

      result.model=-1;
      if (!spectrum["MODSEL"].empty())
      {
        if (spectrum["MODSEL"][1].is_data_line()) result.model=SLHAea::to<int>(spectrum["MODSEL"][1][1]);
        if (spectrum["MODSEL"][3].is_data_line()) result.NMSSM=SLHAea::to<int>(spectrum["MODSEL"][3][1]);
        if (spectrum["MODSEL"][4].is_data_line()) result.RV=SLHAea::to<int>(spectrum["MODSEL"][4][1]);
        if (spectrum["MODSEL"][5].is_data_line()) result.CPV=SLHAea::to<int>(spectrum["MODSEL"][5][1]);
        if (spectrum["MODSEL"][6].is_data_line()) result.FV=SLHAea::to<int>(spectrum["MODSEL"][6][1]);
        if (spectrum["MODSEL"][12].is_data_line()) result.Q=SLHAea::to<double>(spectrum["MODSEL"][12][1]);
      }

      if (result.NMSSM != 0) result.model=result.NMSSM;
      if (result.RV != 0) result.model=-2;
      if (result.CPV != 0) result.model=-2;

      if (!spectrum["SMINPUTS"].empty())
      {
        if (spectrum["SMINPUTS"][1].is_data_line()) result.inv_alpha_em=SLHAea::to<double>(spectrum["SMINPUTS"][1][1]);
        if (spectrum["SMINPUTS"][2].is_data_line()) result.Gfermi=SLHAea::to<double>(spectrum["SMINPUTS"][2][1]);
        if (spectrum["SMINPUTS"][3].is_data_line()) result.alphas_MZ=SLHAea::to<double>(spectrum["SMINPUTS"][3][1]);
        if (spectrum["SMINPUTS"][4].is_data_line()) result.mass_Z=SLHAea::to<double>(spectrum["SMINPUTS"][4][1]);
        if (spectrum["SMINPUTS"][5].is_data_line()) result.mass_b=SLHAea::to<double>(spectrum["SMINPUTS"][5][1]);
        if (spectrum["SMINPUTS"][6].is_data_line()) result.mass_top_pole=SLHAea::to<double>(spectrum["SMINPUTS"][6][1]);
        if (spectrum["SMINPUTS"][7].is_data_line()) result.mass_tau_pole=SLHAea::to<double>(spectrum["SMINPUTS"][7][1]);
        if (spectrum["SMINPUTS"][8].is_data_line()) result.mass_nut=SLHAea::to<double>(spectrum["SMINPUTS"][8][1]);
        if (spectrum["SMINPUTS"][11].is_data_line()) result.mass_e=SLHAea::to<double>(spectrum["SMINPUTS"][11][1]);
        if (spectrum["SMINPUTS"][12].is_data_line()) result.mass_nue=SLHAea::to<double>(spectrum["SMINPUTS"][12][1]);
        if (spectrum["SMINPUTS"][13].is_data_line()) result.mass_mu=SLHAea::to<double>(spectrum["SMINPUTS"][13][1]);
        if (spectrum["SMINPUTS"][14].is_data_line()) result.mass_num=SLHAea::to<double>(spectrum["SMINPUTS"][14][1]);
        if (spectrum["SMINPUTS"][21].is_data_line()) result.mass_d=SLHAea::to<double>(spectrum["SMINPUTS"][21][1]);
        if (spectrum["SMINPUTS"][22].is_data_line()) result.mass_u=SLHAea::to<double>(spectrum["SMINPUTS"][22][1]);
        if (spectrum["SMINPUTS"][23].is_data_line()) result.mass_s=SLHAea::to<double>(spectrum["SMINPUTS"][23][1]);
        if (spectrum["SMINPUTS"][24].is_data_line()) result.mass_c=SLHAea::to<double>(spectrum["SMINPUTS"][24][1]);result.scheme_c_mass=1;
      }

      if (!spectrum["VCKMIN"].empty())
      {
        if (spectrum["VCKMIN"][1].is_data_line()) result.CKM_lambda=SLHAea::to<double>(spectrum["VCKMIN"][1][1]);
        if (spectrum["VCKMIN"][2].is_data_line()) result.CKM_A=SLHAea::to<double>(spectrum["VCKMIN"][2][1]);
        if (spectrum["VCKMIN"][3].is_data_line()) result.CKM_rhobar=SLHAea::to<double>(spectrum["VCKMIN"][3][1]);
        if (spectrum["VCKMIN"][4].is_data_line()) result.CKM_etabar=SLHAea::to<double>(spectrum["VCKMIN"][4][1]);
      }

      if (!spectrum["UPMNSIN"].empty())
      {
        if (spectrum["UPMNSIN"][1].is_data_line()) result.PMNS_theta12=SLHAea::to<double>(spectrum["UPMNSIN"][1][1]);
        if (spectrum["UPMNSIN"][2].is_data_line()) result.PMNS_theta23=SLHAea::to<double>(spectrum["UPMNSIN"][2][1]);
        if (spectrum["UPMNSIN"][3].is_data_line()) result.PMNS_theta13=SLHAea::to<double>(spectrum["UPMNSIN"][3][1]);
        if (spectrum["UPMNSIN"][4].is_data_line()) result.PMNS_delta13=SLHAea::to<double>(spectrum["UPMNSIN"][4][1]);
        if (spectrum["UPMNSIN"][5].is_data_line()) result.PMNS_alpha1=SLHAea::to<double>(spectrum["UPMNSIN"][5][1]);
        if (spectrum["UPMNSIN"][6].is_data_line()) result.PMNS_alpha2=SLHAea::to<double>(spectrum["UPMNSIN"][6][1]);
      }

      if (!spectrum["MINPAR"].empty())
      {
        switch(result.model)
        {
          case 1:
          {
            if (spectrum["MINPAR"][1].is_data_line()) result.m0=SLHAea::to<double>(spectrum["MINPAR"][1][1]);
            if (spectrum["MINPAR"][2].is_data_line()) result.m12=SLHAea::to<double>(spectrum["MINPAR"][2][1]);
            if (spectrum["MINPAR"][3].is_data_line()) result.tan_beta=SLHAea::to<double>(spectrum["MINPAR"][3][1]);
            if (spectrum["MINPAR"][4].is_data_line()) result.sign_mu=SLHAea::to<double>(spectrum["MINPAR"][4][1]);
            if (spectrum["MINPAR"][5].is_data_line()) result.A0=SLHAea::to<double>(spectrum["MINPAR"][5][1]);
            break;
          }
          case 2:
          {
            if (spectrum["MINPAR"][1].is_data_line()) result.Lambda=SLHAea::to<double>(spectrum["MINPAR"][1][1]);
            if (spectrum["MINPAR"][2].is_data_line()) result.Mmess=SLHAea::to<double>(spectrum["MINPAR"][2][1]);
            if (spectrum["MINPAR"][3].is_data_line()) result.tan_beta=SLHAea::to<double>(spectrum["MINPAR"][3][1]);
            if (spectrum["MINPAR"][4].is_data_line()) result.sign_mu=SLHAea::to<double>(spectrum["MINPAR"][4][1]);
            if (spectrum["MINPAR"][5].is_data_line()) result.N5=SLHAea::to<double>(spectrum["MINPAR"][5][1]);
            if (spectrum["MINPAR"][6].is_data_line()) result.cgrav=SLHAea::to<double>(spectrum["MINPAR"][6][1]);
            break;
          }
          case 3:
          {
            if (spectrum["MINPAR"][1].is_data_line()) result.m32=SLHAea::to<double>(spectrum["MINPAR"][1][1]);
            if (spectrum["MINPAR"][2].is_data_line()) result.m0=SLHAea::to<double>(spectrum["MINPAR"][2][1]);
            if (spectrum["MINPAR"][3].is_data_line()) result.tan_beta=SLHAea::to<double>(spectrum["MINPAR"][3][1]);
            if (spectrum["MINPAR"][4].is_data_line()) result.sign_mu=SLHAea::to<double>(spectrum["MINPAR"][4][1]);
            break;
          }
          default:
          {
            if (spectrum["MINPAR"][3].is_data_line()) result.tan_beta=SLHAea::to<double>(spectrum["MINPAR"][3][1]);
            break;
          }
        }
      }

      if (!spectrum["EXTPAR"].empty())
      {
        if (spectrum["EXTPAR"][0].is_data_line()) result.Min=SLHAea::to<double>(spectrum["EXTPAR"][0][1]);
        if (spectrum["EXTPAR"][1].is_data_line()) result.M1_Min=SLHAea::to<double>(spectrum["EXTPAR"][1][1]);
        if (spectrum["EXTPAR"][2].is_data_line()) result.M2_Min=SLHAea::to<double>(spectrum["EXTPAR"][2][1]);
        if (spectrum["EXTPAR"][3].is_data_line()) result.M3_Min=SLHAea::to<double>(spectrum["EXTPAR"][3][1]);
        if (spectrum["EXTPAR"][11].is_data_line()) result.At_Min=SLHAea::to<double>(spectrum["EXTPAR"][11][1]);
        if (spectrum["EXTPAR"][12].is_data_line()) result.Ab_Min=SLHAea::to<double>(spectrum["EXTPAR"][12][1]);
        if (spectrum["EXTPAR"][13].is_data_line()) result.Atau_Min=SLHAea::to<double>(spectrum["EXTPAR"][13][1]);
        if (spectrum["EXTPAR"][21].is_data_line()) result.M2H1_Min=SLHAea::to<double>(spectrum["EXTPAR"][21][1]);
        if (spectrum["EXTPAR"][22].is_data_line()) result.M2H2_Min=SLHAea::to<double>(spectrum["EXTPAR"][22][1]);
        if (spectrum["EXTPAR"][23].is_data_line()) result.mu_Min=SLHAea::to<double>(spectrum["EXTPAR"][23][1]);
        if (spectrum["EXTPAR"][24].is_data_line()) result.M2A_Min=SLHAea::to<double>(spectrum["EXTPAR"][24][1]);
        if (spectrum["EXTPAR"][25].is_data_line()) result.tb_Min=SLHAea::to<double>(spectrum["EXTPAR"][25][1]);
        if (spectrum["EXTPAR"][26].is_data_line()) result.mA_Min=SLHAea::to<double>(spectrum["EXTPAR"][26][1]);
        if (spectrum["EXTPAR"][31].is_data_line()) result.MeL_Min=SLHAea::to<double>(spectrum["EXTPAR"][31][1]);
        if (spectrum["EXTPAR"][32].is_data_line()) result.MmuL_Min=SLHAea::to<double>(spectrum["EXTPAR"][32][1]);
        if (spectrum["EXTPAR"][33].is_data_line()) result.MtauL_Min=SLHAea::to<double>(spectrum["EXTPAR"][33][1]);
        if (spectrum["EXTPAR"][34].is_data_line()) result.MeR_Min=SLHAea::to<double>(spectrum["EXTPAR"][34][1]);
        if (spectrum["EXTPAR"][35].is_data_line()) result.MmuR_Min=SLHAea::to<double>(spectrum["EXTPAR"][35][1]);
        if (spectrum["EXTPAR"][36].is_data_line()) result.MtauR_Min=SLHAea::to<double>(spectrum["EXTPAR"][36][1]);
        if (spectrum["EXTPAR"][41].is_data_line()) result.MqL1_Min=SLHAea::to<double>(spectrum["EXTPAR"][41][1]);
        if (spectrum["EXTPAR"][42].is_data_line()) result.MqL2_Min=SLHAea::to<double>(spectrum["EXTPAR"][42][1]);
        if (spectrum["EXTPAR"][43].is_data_line()) result.MqL3_Min=SLHAea::to<double>(spectrum["EXTPAR"][43][1]);
        if (spectrum["EXTPAR"][44].is_data_line()) result.MuR_Min=SLHAea::to<double>(spectrum["EXTPAR"][44][1]);
        if (spectrum["EXTPAR"][45].is_data_line()) result.McR_Min=SLHAea::to<double>(spectrum["EXTPAR"][45][1]);
        if (spectrum["EXTPAR"][46].is_data_line()) result.MtR_Min=SLHAea::to<double>(spectrum["EXTPAR"][46][1]);
        if (spectrum["EXTPAR"][47].is_data_line()) result.MdR_Min=SLHAea::to<double>(spectrum["EXTPAR"][47][1]);
        if (spectrum["EXTPAR"][48].is_data_line()) result.MsR_Min=SLHAea::to<double>(spectrum["EXTPAR"][48][1]);
        if (spectrum["EXTPAR"][49].is_data_line()) result.MbR_Min=SLHAea::to<double>(spectrum["EXTPAR"][49][1]);
        if (spectrum["EXTPAR"][51].is_data_line()) result.N51=SLHAea::to<double>(spectrum["EXTPAR"][51][1]);
        if (spectrum["EXTPAR"][52].is_data_line()) result.N52=SLHAea::to<double>(spectrum["EXTPAR"][52][1]);
        if (spectrum["EXTPAR"][53].is_data_line()) result.N53=SLHAea::to<double>(spectrum["EXTPAR"][53][1]);
        if (spectrum["EXTPAR"][61].is_data_line()) result.lambdaNMSSM_Min=SLHAea::to<double>(spectrum["EXTPAR"][61][1]);
        if (spectrum["EXTPAR"][62].is_data_line()) result.kappaNMSSM_Min=SLHAea::to<double>(spectrum["EXTPAR"][62][1]);
        if (spectrum["EXTPAR"][63].is_data_line()) result.AlambdaNMSSM_Min=SLHAea::to<double>(spectrum["EXTPAR"][63][1]);
        if (spectrum["EXTPAR"][64].is_data_line()) result.AkappaNMSSM_Min=SLHAea::to<double>(spectrum["EXTPAR"][64][1]);
        if (spectrum["EXTPAR"][65].is_data_line()) result.lambdaSNMSSM_Min=SLHAea::to<double>(spectrum["EXTPAR"][65][1]);
        if (spectrum["EXTPAR"][66].is_data_line()) result.xiFNMSSM_Min=SLHAea::to<double>(spectrum["EXTPAR"][66][1]);
        if (spectrum["EXTPAR"][67].is_data_line()) result.xiSNMSSM_Min=SLHAea::to<double>(spectrum["EXTPAR"][67][1]);
        if (spectrum["EXTPAR"][68].is_data_line()) result.mupNMSSM_Min=SLHAea::to<double>(spectrum["EXTPAR"][68][1]);
        if (spectrum["EXTPAR"][69].is_data_line()) result.mSp2NMSSM_Min=SLHAea::to<double>(spectrum["EXTPAR"][69][1]);
        if (spectrum["EXTPAR"][70].is_data_line()) result.mS2NMSSM_Min=SLHAea::to<double>(spectrum["EXTPAR"][70][1]);
      }

      if (!spectrum["MASS"].empty())
      {
        if (spectrum["MASS"][1].is_data_line()) result.mass_d=SLHAea::to<double>(spectrum["MASS"][1][1]);
        if (spectrum["MASS"][2].is_data_line()) result.mass_u=SLHAea::to<double>(spectrum["MASS"][2][1]);
        if (spectrum["MASS"][3].is_data_line()) result.mass_s=SLHAea::to<double>(spectrum["MASS"][3][1]);
        if (spectrum["MASS"][4].is_data_line()) result.mass_c_pole=SLHAea::to<double>(spectrum["MASS"][4][1]);
        if (spectrum["MASS"][6].is_data_line()) result.mass_t=SLHAea::to<double>(spectrum["MASS"][6][1]);
        if (spectrum["MASS"][11].is_data_line()) result.mass_e=SLHAea::to<double>(spectrum["MASS"][11][1]);
        if (spectrum["MASS"][12].is_data_line()) result.mass_nue=SLHAea::to<double>(spectrum["MASS"][12][1]);
        if (spectrum["MASS"][13].is_data_line()) result.mass_mu=SLHAea::to<double>(spectrum["MASS"][13][1]);
        if (spectrum["MASS"][14].is_data_line()) result.mass_num=SLHAea::to<double>(spectrum["MASS"][14][1]);
        if (spectrum["MASS"][15].is_data_line()) result.mass_tau=result.mass_tau_pole=SLHAea::to<double>(spectrum["MASS"][15][1]);
        if (spectrum["MASS"][16].is_data_line()) result.mass_nut=SLHAea::to<double>(spectrum["MASS"][16][1]);
        if (spectrum["MASS"][21].is_data_line()) result.mass_gluon=SLHAea::to<double>(spectrum["MASS"][21][1]);
        if (spectrum["MASS"][22].is_data_line()) result.mass_photon=SLHAea::to<double>(spectrum["MASS"][22][1]);
        if (spectrum["MASS"][23].is_data_line()) result.mass_Z0=SLHAea::to<double>(spectrum["MASS"][23][1]);
        if (spectrum["MASS"][24].is_data_line()) result.mass_W=SLHAea::to<double>(spectrum["MASS"][24][1]);
        if (spectrum["MASS"][25].is_data_line()) result.mass_h0=SLHAea::to<double>(spectrum["MASS"][25][1]);
        if (spectrum["MASS"][35].is_data_line()) result.mass_H0=SLHAea::to<double>(spectrum["MASS"][35][1]);
        if (spectrum["MASS"][36].is_data_line()) result.mass_A0=SLHAea::to<double>(spectrum["MASS"][36][1]);
        if (spectrum["MASS"][37].is_data_line()) result.mass_H=SLHAea::to<double>(spectrum["MASS"][37][1]);
        if (spectrum["MASS"][39].is_data_line()) result.mass_graviton=SLHAea::to<double>(spectrum["MASS"][39][1]);
        if (spectrum["MASS"][45].is_data_line()) result.mass_H03=SLHAea::to<double>(spectrum["MASS"][45][1]);
        if (spectrum["MASS"][46].is_data_line()) result.mass_A02=SLHAea::to<double>(spectrum["MASS"][46][1]);
        if (spectrum["MASS"][1000001].is_data_line()) result.mass_dnl=SLHAea::to<double>(spectrum["MASS"][1000001][1]);
        if (spectrum["MASS"][1000002].is_data_line()) result.mass_upl=SLHAea::to<double>(spectrum["MASS"][1000002][1]);
        if (spectrum["MASS"][1000003].is_data_line()) result.mass_stl=SLHAea::to<double>(spectrum["MASS"][1000003][1]);
        if (spectrum["MASS"][1000004].is_data_line()) result.mass_chl=SLHAea::to<double>(spectrum["MASS"][1000004][1]);
        if (spectrum["MASS"][1000005].is_data_line()) result.mass_b1=SLHAea::to<double>(spectrum["MASS"][1000005][1]);
        if (spectrum["MASS"][1000006].is_data_line()) result.mass_t1=SLHAea::to<double>(spectrum["MASS"][1000006][1]);
        if (spectrum["MASS"][1000011].is_data_line()) result.mass_el=SLHAea::to<double>(spectrum["MASS"][1000011][1]);
        if (spectrum["MASS"][1000012].is_data_line()) result.mass_nuel=SLHAea::to<double>(spectrum["MASS"][1000012][1]);
        if (spectrum["MASS"][1000013].is_data_line()) result.mass_mul=SLHAea::to<double>(spectrum["MASS"][1000013][1]);
        if (spectrum["MASS"][1000014].is_data_line()) result.mass_numl=SLHAea::to<double>(spectrum["MASS"][1000014][1]);
        if (spectrum["MASS"][1000015].is_data_line()) result.mass_tau1=SLHAea::to<double>(spectrum["MASS"][1000015][1]);
        if (spectrum["MASS"][1000016].is_data_line()) result.mass_nutl=SLHAea::to<double>(spectrum["MASS"][1000016][1]);
        if (spectrum["MASS"][1000021].is_data_line()) result.mass_gluino=SLHAea::to<double>(spectrum["MASS"][1000021][1]);
        if (spectrum["MASS"][1000022].is_data_line()) result.mass_neut[1]=SLHAea::to<double>(spectrum["MASS"][1000022][1]);
        if (spectrum["MASS"][1000023].is_data_line()) result.mass_neut[2]=SLHAea::to<double>(spectrum["MASS"][1000023][1]);
        if (spectrum["MASS"][1000024].is_data_line()) result.mass_cha1=SLHAea::to<double>(spectrum["MASS"][1000024][1]);
        if (spectrum["MASS"][1000025].is_data_line()) result.mass_neut[3]=SLHAea::to<double>(spectrum["MASS"][1000025][1]);
        if (spectrum["MASS"][1000035].is_data_line()) result.mass_neut[4]=SLHAea::to<double>(spectrum["MASS"][1000035][1]);
        if (spectrum["MASS"][1000037].is_data_line()) result.mass_cha2=SLHAea::to<double>(spectrum["MASS"][1000037][1]);
        if (spectrum["MASS"][1000039].is_data_line()) result.mass_gravitino=SLHAea::to<double>(spectrum["MASS"][1000039][1]);
        if (spectrum["MASS"][1000045].is_data_line()) result.mass_neut[5]=SLHAea::to<double>(spectrum["MASS"][1000045][1]);
        if (spectrum["MASS"][2000001].is_data_line()) result.mass_dnr=SLHAea::to<double>(spectrum["MASS"][2000001][1]);
        if (spectrum["MASS"][2000002].is_data_line()) result.mass_upr=SLHAea::to<double>(spectrum["MASS"][2000002][1]);
        if (spectrum["MASS"][2000003].is_data_line()) result.mass_str=SLHAea::to<double>(spectrum["MASS"][2000003][1]);
        if (spectrum["MASS"][2000004].is_data_line()) result.mass_chr=SLHAea::to<double>(spectrum["MASS"][2000004][1]);
        if (spectrum["MASS"][2000005].is_data_line()) result.mass_b2=SLHAea::to<double>(spectrum["MASS"][2000005][1]);
        if (spectrum["MASS"][2000006].is_data_line()) result.mass_t2=SLHAea::to<double>(spectrum["MASS"][2000006][1]);
        if (spectrum["MASS"][2000011].is_data_line()) result.mass_er=SLHAea::to<double>(spectrum["MASS"][2000011][1]);
        if (spectrum["MASS"][2000012].is_data_line()) result.mass_nuer=SLHAea::to<double>(spectrum["MASS"][2000012][1]);
        if (spectrum["MASS"][2000013].is_data_line()) result.mass_mur=SLHAea::to<double>(spectrum["MASS"][2000013][1]);
        if (spectrum["MASS"][2000014].is_data_line()) result.mass_numr=SLHAea::to<double>(spectrum["MASS"][2000014][1]);
        if (spectrum["MASS"][2000015].is_data_line()) result.mass_tau2=SLHAea::to<double>(spectrum["MASS"][2000015][1]);
        if (spectrum["MASS"][2000016].is_data_line()) result.mass_nutr=SLHAea::to<double>(spectrum["MASS"][2000016][1]);
      }

      // The following blocks will only appear for SUSY models so let's not waste time checking them if we're not scanning one of those.
      if (ModelInUse("MSSM63atMGUT") or ModelInUse("MSSM63atQ"))
      {
        // The scale doesn't come through in MODSEL with all spectrum generators
        result.Q = Dep::MSSM_spectrum->get_HE().GetScale();

        if (!spectrum["ALPHA"].empty()) if (spectrum["ALPHA"].back().is_data_line()) result.alpha=SLHAea::to<double>(spectrum["ALPHA"].back().at(0));

        if (!spectrum["STOPMIX"].empty()) for (ie=1;ie<=2;ie++) for (je=1;je<=2;je++)
         if (spectrum["STOPMIX"][max(ie,je)].is_data_line()) result.stop_mix[ie][je]=SLHAea::to<double>(spectrum["STOPMIX"].at(ie,je)[2]);
        if (!spectrum["SBOTMIX"].empty()) for (ie=1;ie<=2;ie++) for (je=1;je<=2;je++)
         if (spectrum["SBOTMIX"][max(ie,je)].is_data_line()) result.sbot_mix[ie][je]=SLHAea::to<double>(spectrum["SBOTMIX"].at(ie,je)[2]);
        if (!spectrum["STAUMIX"].empty()) for (ie=1;ie<=2;ie++) for (je=1;je<=2;je++)
         if (spectrum["STAUMIX"][max(ie,je)].is_data_line()) result.stau_mix[ie][je]=SLHAea::to<double>(spectrum["STAUMIX"].at(ie,je)[2]);
        if (!spectrum["NMIX"].empty()) for (ie=1;ie<=4;ie++) for (je=1;je<=4;je++)
         if (spectrum["NMIX"][max(ie,je)].is_data_line()) result.neut_mix[ie][je]=SLHAea::to<double>(spectrum["NMIX"].at(ie,je)[2]);
        if (!spectrum["NMNMIX"].empty()) for (ie=1;ie<=5;ie++) for (je=1;je<=5;je++)
         if (spectrum["NMNMIX"][max(ie,je)].is_data_line()) result.neut_mix[ie][je]=SLHAea::to<double>(spectrum["NMNMIX"].at(ie,je)[2]);
        if (!spectrum["UMIX"].empty()) for (ie=1;ie<=2;ie++) for (je=1;je<=2;je++)
         if (spectrum["UMIX"][max(ie,je)].is_data_line()) result.charg_Umix[ie][je]=SLHAea::to<double>(spectrum["UMIX"].at(ie,je)[2]);
        if (!spectrum["VMIX"].empty()) for (ie=1;ie<=2;ie++) for (je=1;je<=2;je++)
         if (spectrum["VMIX"][max(ie,je)].is_data_line()) result.charg_Vmix[ie][je]=SLHAea::to<double>(spectrum["VMIX"].at(ie,je)[2]);

        if (!spectrum["GAUGE"].empty())
        {
          if (spectrum["GAUGE"][1].is_data_line()) result.gp_Q=SLHAea::to<double>(spectrum["GAUGE"][1][1]);
          if (spectrum["GAUGE"][2].is_data_line()) result.g2_Q=SLHAea::to<double>(spectrum["GAUGE"][2][1]);
          if (spectrum["GAUGE"][3].is_data_line()) result.g3_Q=SLHAea::to<double>(spectrum["GAUGE"][3][1]);
        }

        if (!spectrum["YU"].empty()) for (ie=1;ie<=3;ie++) if (spectrum["YU"][ie].is_data_line()) result.yut[ie]=SLHAea::to<double>(spectrum["YU"].at(ie,ie)[2]);
        if (!spectrum["YD"].empty()) for (ie=1;ie<=3;ie++) if (spectrum["YD"][ie].is_data_line()) result.yub[ie]=SLHAea::to<double>(spectrum["YD"].at(ie,ie)[2]);
        if (!spectrum["YE"].empty()) for (ie=1;ie<=3;ie++) if (spectrum["YE"][ie].is_data_line()) result.yutau[ie]=SLHAea::to<double>(spectrum["YE"].at(ie,ie)[2]);

        if (!spectrum["HMIX"].empty())
        {
          if (spectrum["HMIX"][1].is_data_line()) result.mu_Q=SLHAea::to<double>(spectrum["HMIX"][1][1]);
          if (spectrum["HMIX"][2].is_data_line()) result.tanb_GUT=SLHAea::to<double>(spectrum["HMIX"][2][1]);
          if (spectrum["HMIX"][3].is_data_line()) result.Higgs_VEV=SLHAea::to<double>(spectrum["HMIX"][3][1]);
          if (spectrum["HMIX"][4].is_data_line()) result.mA2_Q=SLHAea::to<double>(spectrum["HMIX"][4][1]);
        }

        if (!spectrum["NMHMIX"].empty()) for (ie=1;ie<=3;ie++) for (je=1;je<=3;je++)
         if (spectrum["NMHMIX"][max(ie,je)].is_data_line()) result.H0_mix[ie][je]=SLHAea::to<double>(spectrum["NMHMIX"].at(ie,je)[2]);

        if (!spectrum["NMAMIX"].empty()) for (ie=1;ie<=2;ie++) for (je=1;je<=2;je++)
         if (spectrum["NMAMIX"][max(ie,je)].is_data_line()) result.A0_mix[ie][je]=SLHAea::to<double>(spectrum["NMAMIX"].at(ie,je)[2]);

        if (!spectrum["MSOFT"].empty())
        {
          if (!spectrum["MSOFT"].front().empty()) result.MSOFT_Q=SLHAea::to<double>(spectrum["MSOFT"].front().at(3));
          if (spectrum["MSOFT"][1].is_data_line()) result.M1_Q=SLHAea::to<double>(spectrum["MSOFT"][1][1]);
          if (spectrum["MSOFT"][2].is_data_line()) result.M2_Q=SLHAea::to<double>(spectrum["MSOFT"][2][1]);
          if (spectrum["MSOFT"][3].is_data_line()) result.M3_Q=SLHAea::to<double>(spectrum["MSOFT"][3][1]);
          if (spectrum["MSOFT"][21].is_data_line()) result.M2H1_Q=SLHAea::to<double>(spectrum["MSOFT"][21][1]);
          if (spectrum["MSOFT"][22].is_data_line()) result.M2H2_Q=SLHAea::to<double>(spectrum["MSOFT"][22][1]);
          if (spectrum["MSOFT"][31].is_data_line()) result.MeL_Q=SLHAea::to<double>(spectrum["MSOFT"][31][1]);
          if (spectrum["MSOFT"][32].is_data_line()) result.MmuL_Q=SLHAea::to<double>(spectrum["MSOFT"][32][1]);
          if (spectrum["MSOFT"][33].is_data_line()) result.MtauL_Q=SLHAea::to<double>(spectrum["MSOFT"][33][1]);
          if (spectrum["MSOFT"][34].is_data_line()) result.MeR_Q=SLHAea::to<double>(spectrum["MSOFT"][34][1]);
          if (spectrum["MSOFT"][35].is_data_line()) result.MmuR_Q=SLHAea::to<double>(spectrum["MSOFT"][35][1]);
          if (spectrum["MSOFT"][36].is_data_line()) result.MtauR_Q=SLHAea::to<double>(spectrum["MSOFT"][36][1]);
          if (spectrum["MSOFT"][41].is_data_line()) result.MqL1_Q=SLHAea::to<double>(spectrum["MSOFT"][41][1]);
          if (spectrum["MSOFT"][42].is_data_line()) result.MqL2_Q=SLHAea::to<double>(spectrum["MSOFT"][42][1]);
          if (spectrum["MSOFT"][43].is_data_line()) result.MqL3_Q=SLHAea::to<double>(spectrum["MSOFT"][43][1]);
          if (spectrum["MSOFT"][44].is_data_line()) result.MuR_Q=SLHAea::to<double>(spectrum["MSOFT"][44][1]);
          if (spectrum["MSOFT"][45].is_data_line()) result.McR_Q=SLHAea::to<double>(spectrum["MSOFT"][45][1]);
          if (spectrum["MSOFT"][46].is_data_line()) result.MtR_Q=SLHAea::to<double>(spectrum["MSOFT"][46][1]);
          if (spectrum["MSOFT"][47].is_data_line()) result.MdR_Q=SLHAea::to<double>(spectrum["MSOFT"][47][1]);
          if (spectrum["MSOFT"][48].is_data_line()) result.MsR_Q=SLHAea::to<double>(spectrum["MSOFT"][48][1]);
          if (spectrum["MSOFT"][49].is_data_line()) result.MbR_Q=SLHAea::to<double>(spectrum["MSOFT"][49][1]);
        }

        if (!spectrum["AU"].empty())
        {
          if (spectrum["AU"][1].is_data_line()) result.A_u=SLHAea::to<double>(spectrum["AU"].at(1,1)[2]);
          if (spectrum["AU"][2].is_data_line()) result.A_c=SLHAea::to<double>(spectrum["AU"].at(2,2)[2]);
          if (spectrum["AU"][3].is_data_line()) result.A_t=SLHAea::to<double>(spectrum["AU"].at(3,3)[2]);
        }

        if (!spectrum["AD"].empty())
        {
          if (spectrum["AD"][1].is_data_line()) result.A_d=SLHAea::to<double>(spectrum["AD"].at(1,1)[2]);
          if (spectrum["AD"][2].is_data_line()) result.A_s=SLHAea::to<double>(spectrum["AD"].at(2,2)[2]);
          if (spectrum["AD"][3].is_data_line()) result.A_b=SLHAea::to<double>(spectrum["AD"].at(3,3)[2]);
        }

        if (!spectrum["AE"].empty())
        {
          if (spectrum["AE"][1].is_data_line()) result.A_e=SLHAea::to<double>(spectrum["AE"].at(1,1)[2]);
          if (spectrum["AE"][2].is_data_line()) result.A_mu=SLHAea::to<double>(spectrum["AE"].at(2,2)[2]);
          if (spectrum["AE"][3].is_data_line()) result.A_tau=SLHAea::to<double>(spectrum["AE"].at(3,3)[2]);
        }

        if (!spectrum["NMSSMRUN"].empty())
        {
          if (spectrum["NMSSMRUN"][1].is_data_line()) result.lambdaNMSSM=SLHAea::to<double>(spectrum["NMSSMRUN"][1][1]);
          if (spectrum["NMSSMRUN"][2].is_data_line()) result.kappaNMSSM=SLHAea::to<double>(spectrum["NMSSMRUN"][2][1]);
          if (spectrum["NMSSMRUN"][3].is_data_line()) result.AlambdaNMSSM=SLHAea::to<double>(spectrum["NMSSMRUN"][3][1]);
          if (spectrum["NMSSMRUN"][4].is_data_line()) result.AkappaNMSSM=SLHAea::to<double>(spectrum["NMSSMRUN"][4][1]);
          if (spectrum["NMSSMRUN"][5].is_data_line()) result.lambdaSNMSSM=SLHAea::to<double>(spectrum["NMSSMRUN"][5][1]);
          if (spectrum["NMSSMRUN"][6].is_data_line()) result.xiFNMSSM=SLHAea::to<double>(spectrum["NMSSMRUN"][6][1]);
          if (spectrum["NMSSMRUN"][7].is_data_line()) result.xiSNMSSM=SLHAea::to<double>(spectrum["NMSSMRUN"][7][1]);
          if (spectrum["NMSSMRUN"][8].is_data_line()) result.mupNMSSM=SLHAea::to<double>(spectrum["NMSSMRUN"][8][1]);
          if (spectrum["NMSSMRUN"][9].is_data_line()) result.mSp2NMSSM=SLHAea::to<double>(spectrum["NMSSMRUN"][9][1]);
          if (spectrum["NMSSMRUN"][10].is_data_line()) result.mS2NMSSM=SLHAea::to<double>(spectrum["NMSSMRUN"][10][1]);
        }

        if (!spectrum["USQMIX"].empty()) for (ie=1;ie<=6;ie++) for (je=1;je<=6;je++)
         if (spectrum["USQMIX"][max(ie,je)].is_data_line()) result.sU_mix[ie][je]=SLHAea::to<double>(spectrum["USQMIX"].at(ie,je)[2]);
        if (!spectrum["DSQMIX"].empty()) for (ie=1;ie<=6;ie++) for (je=1;je<=6;je++)
         if (spectrum["DSQMIX"][max(ie,je)].is_data_line()) result.sD_mix[ie][je]=SLHAea::to<double>(spectrum["DSQMIX"].at(ie,je)[2]);
        if (!spectrum["SELMIX"].empty()) for (ie=1;ie<=6;ie++) for (je=1;je<=6;je++)
         if (spectrum["SELMIX"][max(ie,je)].is_data_line()) result.sE_mix[ie][je]=SLHAea::to<double>(spectrum["SELMIX"].at(ie,je)[2]);
        if (!spectrum["SNUMIX"].empty()) for (ie=1;ie<=3;ie++) for (je=1;je<=3;je++)
         if (spectrum["SNUMIX"][max(ie,je)].is_data_line()) result.sNU_mix[ie][je]=SLHAea::to<double>(spectrum["SNUMIX"].at(ie,je)[2]);

        if (!spectrum["MSQ2"].empty()) for (ie=1;ie<=3;ie++) for (je=1;je<=3;je++)
         if (spectrum["MSQ2"][max(ie,je)].is_data_line()) result.sCKM_msq2[ie][je]=SLHAea::to<double>(spectrum["MSQ2"].at(ie,je)[2]);
        if (!spectrum["MSL2"].empty()) for (ie=1;ie<=3;ie++) for (je=1;je<=3;je++)
         if (spectrum["MSL2"][max(ie,je)].is_data_line()) result.sCKM_msl2[ie][je]=SLHAea::to<double>(spectrum["MSL2"].at(ie,je)[2]);
        if (!spectrum["MSD2"].empty()) for (ie=1;ie<=3;ie++) for (je=1;je<=3;je++)
         if (spectrum["MSD2"][max(ie,je)].is_data_line()) result.sCKM_msd2[ie][je]=SLHAea::to<double>(spectrum["MSD2"].at(ie,je)[2]);
        if (!spectrum["MSU2"].empty()) for (ie=1;ie<=3;ie++) for (je=1;je<=3;je++)
         if (spectrum["MSU2"][max(ie,je)].is_data_line()) result.sCKM_msu2[ie][je]=SLHAea::to<double>(spectrum["MSU2"].at(ie,je)[2]);
        if (!spectrum["MSE2"].empty()) for (ie=1;ie<=3;ie++) for (je=1;je<=3;je++)
         if (spectrum["MSE2"][max(ie,je)].is_data_line()) result.sCKM_mse2[ie][je]=SLHAea::to<double>(spectrum["MSE2"].at(ie,je)[2]);

        if (!spectrum["IMVCKM"].empty()) for (ie=1;ie<=3;ie++) for (je=1;je<=3;je++)
         if (spectrum["IMVCKM"][max(ie,je)].is_data_line()) result.IMCKM[ie][je]=SLHAea::to<double>(spectrum["IMVCKM"].at(ie,je)[2]);
        if (!spectrum["IMVCKM"].empty()) for (ie=1;ie<=3;ie++) for (je=1;je<=3;je++)
         if (spectrum["IMVCKM"][max(ie,je)].is_data_line()) result.IMCKM[ie][je]=SLHAea::to<double>(spectrum["IMVCKM"].at(ie,je)[2]);

        if (!spectrum["UPMNS"].empty()) for (ie=1;ie<=3;ie++) for (je=1;je<=3;je++)
         if (spectrum["UPMNS"][max(ie,je)].is_data_line()) result.PMNS_U[ie][je]=SLHAea::to<double>(spectrum["UPMNS"].at(ie,je)[2]);

        if (!spectrum["TU"].empty()) for (ie=1;ie<=3;ie++) for (je=1;je<=3;je++)
         if (spectrum["TU"][max(ie,je)].is_data_line()) result.TU[ie][je]=SLHAea::to<double>(spectrum["TU"].at(ie,je)[2]);
        if (!spectrum["TD"].empty()) for (ie=1;ie<=3;ie++) for (je=1;je<=3;je++)
         if (spectrum["TD"][max(ie,je)].is_data_line()) result.TD[ie][je]=SLHAea::to<double>(spectrum["TD"].at(ie,je)[2]);
        if (!spectrum["TE"].empty()) for (ie=1;ie<=3;ie++) for (je=1;je<=3;je++)
         if (spectrum["TE"][max(ie,je)].is_data_line()) result.TE[ie][je]=SLHAea::to<double>(spectrum["TE"].at(ie,je)[2]);
      }

      else if (ModelInUse("WC"))
      {
        // The Higgs mass doesn't come through in the SLHAea object, as that's only for SLHA2 SM inputs.
        result.mass_h0 = Dep::SM_spectrum->get(Par::Pole_Mass, "h0_1");
        // Set the scale.
        result.Q = result.mass_Z;
      }

	  if(byVal(result.mass_c_pole)>0.&&byVal(result.scheme_c_mass)<0)
	  {
		if(byVal(result.mass_c_pole)<1.5) result.mass_c=BEreq::mcmc_from_pole(byVal(result.mass_c_pole),1,&result);
		else if(byVal(result.mass_c_pole)<1.75) result.mass_c=BEreq::mcmc_from_pole(byVal(result.mass_c_pole),2,&result);
		else result.mass_c=BEreq::mcmc_from_pole(byVal(result.mass_c_pole),3,&result);
	  }

      BEreq::slha_adjust(&result);

      // Set the Z and W widths
      result.width_Z = Dep::Z_decay_rates->width_in_GeV;
      result.width_W = Dep::W_plus_decay_rates->width_in_GeV;

      // If requested, override the SuperIso b pole mass with the SpecBit value and recompute the 1S b mass.
      if (runOptions->getValueOrDef<bool>(false, "take_b_pole_mass_from_spectrum"))
      {
        if (ModelInUse("MSSM63atMGUT") or ModelInUse("MSSM63atQ"))
        {
          result.mass_h0 = Dep::MSSM_spectrum->get(Par::Pole_Mass, "h0_1");
        }
        else if (ModelInUse("WC"))
        {
          result.mass_h0 = Dep::SM_spectrum->get(Par::Pole_Mass, "h0_1");
        }
        result.mass_b_1S = BEreq::mb_1S(&result);
      }

      if (ModelInUse("WC"))
      {

        // Tell SuperIso to do its Wilson coefficient calculations for the SM.
        // We will adjust them with our BSM deviations in backend convenience
        // functions before we send them to SuperIso's observable calculation functions.
        result.SM = 1;

        // So far our model only deals with 5 operators: O_7, O_9, O_10, Q_1 and Q_2.
        result.Re_DeltaC7  = *Param["Re_DeltaC7"];
        result.Im_DeltaC7  = *Param["Im_DeltaC7"];
        result.Re_DeltaC9  = *Param["Re_DeltaC9"];
        result.Im_DeltaC9  = *Param["Im_DeltaC9"];
        result.Re_DeltaC10 = *Param["Re_DeltaC10"];
        result.Im_DeltaC10 = *Param["Im_DeltaC10"];
        result.Re_DeltaCQ1 = *Param["Re_DeltaCQ1"];
        result.Im_DeltaCQ1 = *Param["Im_DeltaCQ1"];
        result.Re_DeltaCQ2 = *Param["Re_DeltaCQ2"];
        result.Im_DeltaCQ2 = *Param["Im_DeltaCQ2"];

        /* Lines below are valid only in the flavour universal case
           deltaC[1..10] = Cmu[1..10], deltaC[11..20] = Ce[1..10], deltaC[21..30] = Ctau[1..10]
           deltaCQ[1,2] = CQmu[1,2], deltaCQ[1,2] = CQe[1,2], deltaCQ[1,2] = CQtau[1,2] */

        result.deltaC[7]=result.deltaC[17]=result.deltaC[27]=std::complex<double>(result.Re_DeltaC7, result.Im_DeltaC7);
        result.deltaC[9]=result.deltaC[19]=result.deltaC[29]=std::complex<double>(result.Re_DeltaC9, result.Im_DeltaC9);
        result.deltaC[10]=result.deltaC[20]=result.deltaC[30]=std::complex<double>(result.Re_DeltaC10, result.Im_DeltaC10);

        result.deltaCQ[1]=result.deltaCQ[3]=result.deltaCQ[5]=std::complex<double>(result.Re_DeltaCQ1, result.Im_DeltaCQ1);
        result.deltaCQ[2]=result.deltaCQ[4]=result.deltaCQ[6]=std::complex<double>(result.Re_DeltaCQ2, result.Im_DeltaCQ2);
      }

      if (flav_debug) cout<<"Finished SI_fill"<<endl;
    }

    /// NEW! Fill SuperIso nuisance structure
    void SI_nuisance_fill(nuisance &nuislist)
    {
      using namespace Pipes::SI_nuisance_fill;
	  if (flav_debug) cout<<"Starting SI_nuisance_fill"<<endl;

      parameters const& param = *Dep::SuperIso_modelinfo;

	  BEreq::set_nuisance(&nuislist);
	  BEreq::set_nuisance_value_from_param(&nuislist,&param);

	  /* Here the nuisance parameters which should not be used for the correlation calculation have to be given a zero standard deviation.
	     E.g. nuislist.mass_b.dev=0.; */

      if (flav_debug) cout<<"Finished SI_nuisance_fill"<<endl;
    }

    /// NEW! Compute values of list of observables
    void SI_compute_obs_list(SI_observable_map& result)  // TO BE MODIFIED
    {
      using namespace Pipes::SI_compute_obs_list;
      if (flav_debug) cout<<"Starting SI_compute_obs_list"<<endl;

      const parameters& param = *Dep::SuperIso_modelinfo;
      const nuisance& nuislist = *Dep::SuperIso_nuisance;
      const std::vector<std::string>& obslist = runOptions->getValue<std::vector<std::string>>("SuperIso_obs_list");

      // --- Needed for SuperIso backend
      int nObservables = obslist.size();

      char obsnames[nObservables][50];
      for(int iObservable = 0; iObservable < nObservables; iObservable++) {
          strcpy(obsnames[iObservable], obslist[iObservable].c_str());
      }

      double *res;
      // Reserve memory
      res = (double *) calloc(nObservables, sizeof(double));
      // --- Needed for SuperIso backend

      BEreq::get_predictions_nuisance((char**)obsnames, &nObservables, &res, &param, &nuislist);

      for(int iObservable = 0; iObservable < nObservables; ++iObservable) {
          result[obslist[iObservable]] = res[iObservable];
      }

      // Free memory
      free(res);
      if (flav_debug) {
          for(int iObservable = 0; iObservable < nObservables; ++iObservable) {
              printf("%s=%.4e\n", obsnames[iObservable], result[obslist[iObservable]]);
          }
      }

      if (flav_debug) {
          cout<<"Finished SI_compute_obs_list"<<endl;
      }
    }

    /// NEW! Compute covariance matrix for a list of observables
    void SI_theory_covariance(SI_covariance_map& result)  // TO BE MODIFIED
    {
      using namespace Pipes::SI_theory_covariance;
      if (flav_debug) cout<<"Starting SI_theory_covariance"<<endl;

      const parameters& param = *Dep::SuperIso_modelinfo;
      const nuisance& nuislist = *Dep::SuperIso_nuisance;
      const std::vector<std::string>& obslist = runOptions->getValue<std::vector<std::string>>("SuperIso_obs_list");

      // --- Needed for SuperIso backend
      int nObservables = obslist.size();

      char obsnames[nObservables][50];
      for(int iObservable = 0; iObservable < nObservables; ++iObservable) {
          strcpy(obsnames[iObservable], obslist[iObservable].c_str());
      }

      int nNuisance=161;
      char namenuisance[nNuisance+1][50];
      BEreq::observables(0, NULL, 0, NULL, NULL, &nuislist, (char **)namenuisance, &param); // Initialization of namenuisance

      const int ncorrnuis=463;
      nuiscorr corrnuis[ncorrnuis]={ // List of nuisance correlations, below between the form factors
      {"a00_BK","a10_BK",-0.39},
      {"a00_BK","a20_BK",-0.71},
      {"a00_BK","a30_BK",-0.63},
      {"a00_BK","a0p_BK",0.49},
      {"a00_BK","a1p_BK",-0.03},
      {"a00_BK","a2p_BK",-0.22},
      {"a00_BK","a0T_BK",0.16},
      {"a00_BK","a1T_BK",-0.08},
      {"a00_BK","a2T_BK",-0.09},
      {"a10_BK","a20_BK",0.66},
      {"a10_BK","a30_BK",0.26},
      {"a10_BK","a0p_BK",0.05},
      {"a10_BK","a1p_BK",0.72},
      {"a10_BK","a2p_BK",0.48},
      {"a10_BK","a0T_BK",-0.08},
      {"a10_BK","a1T_BK",0.03},
      {"a10_BK","a2T_BK",0.01},
      {"a20_BK","a30_BK",0.54},
      {"a20_BK","a0p_BK",-0.17},
      {"a20_BK","a1p_BK",0.51},
      {"a20_BK","a2p_BK",0.59},
      {"a20_BK","a0T_BK",-0.16},
      {"a20_BK","a1T_BK",0.05},
      {"a20_BK","a2T_BK",0.09},
      {"a30_BK","a0p_BK",0.05},
      {"a30_BK","a1p_BK",0.14},
      {"a30_BK","a2p_BK",0.05},
      {"a30_BK","a1T_BK",0.03},
      {"a30_BK","a2T_BK",-0.01},
      {"a0p_BK","a1p_BK",0.09},
      {"a0p_BK","a2p_BK",-0.47},
      {"a0p_BK","a0T_BK",0.34},
      {"a0p_BK","a1T_BK",-0.06},
      {"a0p_BK","a2T_BK",-0.28},
      {"a1p_BK","a2p_BK",0.43},
      {"a1p_BK","a0T_BK",-0.06},
      {"a1p_BK","a1T_BK",0.11},
      {"a1p_BK","a2T_BK",-0.04},
      {"a2p_BK","a0T_BK",-0.32},
      {"a2p_BK","a1T_BK",-0.05},
      {"a2p_BK","a2T_BK",0.29},
      {"a0T_BK","a2T_BK",-0.35},
      {"a1T_BK","a2T_BK",0.21},
      {"a0A0_BKstar","a1A0_BKstar",0.633689},
      {"a0A0_BKstar","a2A0_BKstar",0.0575305},
      {"a0A0_BKstar","a0A1_BKstar",0.0615373},
      {"a0A0_BKstar","a1A1_BKstar",0.169044},
      {"a0A0_BKstar","a2A1_BKstar",0.240647},
      {"a0A0_BKstar","a0A12_BKstar",0.999},
      {"a0A0_BKstar","a1A12_BKstar",0.609495},
      {"a0A0_BKstar","a2A12_BKstar",0.0355029},
      {"a0A0_BKstar","a0V_BKstar",0.0393316},
      {"a0A0_BKstar","a1V_BKstar",0.178168},
      {"a0A0_BKstar","a2V_BKstar",0.195421},
      {"a0A0_BKstar","a0T1_BKstar",0.0164883},
      {"a0A0_BKstar","a1T1_BKstar",0.0452644},
      {"a0A0_BKstar","a2T1_BKstar",-0.158457},
      {"a0A0_BKstar","a0T2_BKstar",0.0217598},
      {"a0A0_BKstar","a1T2_BKstar",0.150264},
      {"a0A0_BKstar","a2T2_BKstar",0.0311554},
      {"a0A0_BKstar","a0T23_BKstar",0.451075},
      {"a0A0_BKstar","a1T23_BKstar",0.104911},
      {"a0A0_BKstar","a2T23_BKstar",-0.347452},
      {"a1A0_BKstar","a2A0_BKstar",0.486314},
      {"a1A0_BKstar","a0A1_BKstar",0.300568},
      {"a1A0_BKstar","a1A1_BKstar",0.472167},
      {"a1A0_BKstar","a2A1_BKstar",0.321569},
      {"a1A0_BKstar","a0A12_BKstar",0.621073},
      {"a1A0_BKstar","a1A12_BKstar",0.756564},
      {"a1A0_BKstar","a2A12_BKstar",0.46913},
      {"a1A0_BKstar","a0V_BKstar",0.280475},
      {"a1A0_BKstar","a1V_BKstar",0.553833},
      {"a1A0_BKstar","a2V_BKstar",0.0752289},
      {"a1A0_BKstar","a0T1_BKstar",0.241261},
      {"a1A0_BKstar","a1T1_BKstar",0.448732},
      {"a1A0_BKstar","a2T1_BKstar",-0.182566},
      {"a1A0_BKstar","a0T2_BKstar",0.24405},
      {"a1A0_BKstar","a1T2_BKstar",0.495428},
      {"a1A0_BKstar","a2T2_BKstar",0.181197},
      {"a1A0_BKstar","a0T23_BKstar",0.573536},
      {"a1A0_BKstar","a1T23_BKstar",0.53356},
      {"a1A0_BKstar","a2T23_BKstar",-0.19402},
      {"a2A0_BKstar","a0A1_BKstar",0.515475},
      {"a2A0_BKstar","a1A1_BKstar",0.743538},
      {"a2A0_BKstar","a2A1_BKstar",0.661689},
      {"a2A0_BKstar","a0A12_BKstar",0.05503},
      {"a2A0_BKstar","a1A12_BKstar",0.326536},
      {"a2A0_BKstar","a2A12_BKstar",0.665149},
      {"a2A0_BKstar","a0V_BKstar",0.545574},
      {"a2A0_BKstar","a1V_BKstar",0.731129},
      {"a2A0_BKstar","a2V_BKstar",0.00860747},
      {"a2A0_BKstar","a0T1_BKstar",0.475401},
      {"a2A0_BKstar","a1T1_BKstar",0.674459},
      {"a2A0_BKstar","a2T1_BKstar",-0.125613},
      {"a2A0_BKstar","a0T2_BKstar",0.474921},
      {"a2A0_BKstar","a1T2_BKstar",0.701034},
      {"a2A0_BKstar","a2T2_BKstar",0.273928},
      {"a2A0_BKstar","a0T23_BKstar",0.663092},
      {"a2A0_BKstar","a1T23_BKstar",0.426721},
      {"a2A0_BKstar","a2T23_BKstar",-0.318126},
      {"a0A1_BKstar","a1A1_BKstar",0.702494},
      {"a0A1_BKstar","a2A1_BKstar",0.153604},
      {"a0A1_BKstar","a0A12_BKstar",0.0555181},
      {"a0A1_BKstar","a1A12_BKstar",0.0105278},
      {"a0A1_BKstar","a2A12_BKstar",0.205325},
      {"a0A1_BKstar","a0V_BKstar",0.922982},
      {"a0A1_BKstar","a1V_BKstar",0.731071},
      {"a0A1_BKstar","a2V_BKstar",-0.35833},
      {"a0A1_BKstar","a0T1_BKstar",0.899902},
      {"a0A1_BKstar","a1T1_BKstar",0.68101},
      {"a0A1_BKstar","a2T1_BKstar",-0.468536},
      {"a0A1_BKstar","a0T2_BKstar",0.899616},
      {"a0A1_BKstar","a1T2_BKstar",0.692637},
      {"a0A1_BKstar","a2T2_BKstar",-0.202092},
      {"a0A1_BKstar","a0T23_BKstar",0.564524},
      {"a0A1_BKstar","a1T23_BKstar",0.30614},
      {"a0A1_BKstar","a2T23_BKstar",-0.362387},
      {"a1A1_BKstar","a2A1_BKstar",0.747682},
      {"a1A1_BKstar","a0A12_BKstar",0.166407},
      {"a1A1_BKstar","a1A12_BKstar",0.284455},
      {"a1A1_BKstar","a2A12_BKstar",0.480939},
      {"a1A1_BKstar","a0V_BKstar",0.678427},
      {"a1A1_BKstar","a1V_BKstar",0.862475},
      {"a1A1_BKstar","a2V_BKstar",-0.0633373},
      {"a1A1_BKstar","a0T1_BKstar",0.630965},
      {"a1A1_BKstar","a1T1_BKstar",0.760031},
      {"a1A1_BKstar","a2T1_BKstar",-0.373063},
      {"a1A1_BKstar","a0T2_BKstar",0.634299},
      {"a1A1_BKstar","a1T2_BKstar",0.915195},
      {"a1A1_BKstar","a2T2_BKstar",0.226359},
      {"a1A1_BKstar","a0T23_BKstar",0.695868},
      {"a1A1_BKstar","a1T23_BKstar",0.457217},
      {"a1A1_BKstar","a2T23_BKstar",-0.440385},
      {"a2A1_BKstar","a0A12_BKstar",0.245093},
      {"a2A1_BKstar","a1A12_BKstar",0.316222},
      {"a2A1_BKstar","a2A12_BKstar",0.535592},
      {"a2A1_BKstar","a0V_BKstar",0.194572},
      {"a2A1_BKstar","a1V_BKstar",0.517205},
      {"a2A1_BKstar","a2V_BKstar",0.253395},
      {"a2A1_BKstar","a0T1_BKstar",0.135433},
      {"a2A1_BKstar","a1T1_BKstar",0.400555},
      {"a2A1_BKstar","a2T1_BKstar",-0.127614},
      {"a2A1_BKstar","a0T2_BKstar",0.140073},
      {"a2A1_BKstar","a1T2_BKstar",0.59715},
      {"a2A1_BKstar","a2T2_BKstar",0.402339},
      {"a2A1_BKstar","a0T23_BKstar",0.553618},
      {"a2A1_BKstar","a1T23_BKstar",0.252273},
      {"a2A1_BKstar","a2T23_BKstar",-0.40495},
      {"a0A12_BKstar","a1A12_BKstar",0.617726},
      {"a0A12_BKstar","a2A12_BKstar",0.0443495},
      {"a0A12_BKstar","a0V_BKstar",0.0335119},
      {"a0A12_BKstar","a1V_BKstar",0.172545},
      {"a0A12_BKstar","a2V_BKstar",0.202654},
      {"a0A12_BKstar","a0T1_BKstar",0.010524},
      {"a0A12_BKstar","a1T1_BKstar",0.0391092},
      {"a0A12_BKstar","a2T1_BKstar",-0.154503},
      {"a0A12_BKstar","a0T2_BKstar",0.0157345},
      {"a0A12_BKstar","a1T2_BKstar",0.145566},
      {"a0A12_BKstar","a2T2_BKstar",0.0343108},
      {"a0A12_BKstar","a0T23_BKstar",0.448529},
      {"a0A12_BKstar","a1T23_BKstar",0.100428},
      {"a0A12_BKstar","a2T23_BKstar",-0.348189},
      {"a1A12_BKstar","a2A12_BKstar",0.700096},
      {"a1A12_BKstar","a0V_BKstar",-0.0140921},
      {"a1A12_BKstar","a1V_BKstar",0.350603},
      {"a1A12_BKstar","a2V_BKstar",0.234661},
      {"a1A12_BKstar","a0T1_BKstar",-0.0243359},
      {"a1A12_BKstar","a1T1_BKstar",0.234436},
      {"a1A12_BKstar","a2T1_BKstar",-0.0609604},
      {"a1A12_BKstar","a0T2_BKstar",-0.0214341},
      {"a1A12_BKstar","a1T2_BKstar",0.28577},
      {"a1A12_BKstar","a2T2_BKstar",0.248484},
      {"a1A12_BKstar","a0T23_BKstar",0.386552},
      {"a1A12_BKstar","a1T23_BKstar",0.537224},
      {"a1A12_BKstar","a2T23_BKstar",-0.0377692},
      {"a2A12_BKstar","a0V_BKstar",0.196602},
      {"a2A12_BKstar","a1V_BKstar",0.485841},
      {"a2A12_BKstar","a2V_BKstar",0.139392},
      {"a2A12_BKstar","a0T1_BKstar",0.173373},
      {"a2A12_BKstar","a1T1_BKstar",0.402833},
      {"a2A12_BKstar","a2T1_BKstar",-0.0779755},
      {"a2A12_BKstar","a0T2_BKstar",0.173504},
      {"a2A12_BKstar","a1T2_BKstar",0.422339},
      {"a2A12_BKstar","a2T2_BKstar",0.253397},
      {"a2A12_BKstar","a0T23_BKstar",0.403006},
      {"a2A12_BKstar","a1T23_BKstar",0.571118},
      {"a2A12_BKstar","a2T23_BKstar",-0.0276853},
      {"a0V_BKstar","a1V_BKstar",0.757379},
      {"a0V_BKstar","a2V_BKstar",-0.397005},
      {"a0V_BKstar","a0T1_BKstar",0.901557},
      {"a0V_BKstar","a1T1_BKstar",0.716437},
      {"a0V_BKstar","a2T1_BKstar",-0.400738},
      {"a0V_BKstar","a0T2_BKstar",0.899723},
      {"a0V_BKstar","a1T2_BKstar",0.702845},
      {"a0V_BKstar","a2T2_BKstar",-0.147198},
      {"a0V_BKstar","a0T23_BKstar",0.573778},
      {"a0V_BKstar","a1T23_BKstar",0.277866},
      {"a0V_BKstar","a2T23_BKstar",-0.360324},
      {"a1V_BKstar","a2V_BKstar",0.0346143},
      {"a1V_BKstar","a0T1_BKstar",0.709376},
      {"a1V_BKstar","a1T1_BKstar",0.906137},
      {"a1V_BKstar","a2T1_BKstar",-0.236675},
      {"a1V_BKstar","a0T2_BKstar",0.709471},
      {"a1V_BKstar","a1T2_BKstar",0.916209},
      {"a1V_BKstar","a2T2_BKstar",0.261192},
      {"a1V_BKstar","a0T23_BKstar",0.70266},
      {"a1V_BKstar","a1T23_BKstar",0.480664},
      {"a1V_BKstar","a2T23_BKstar",-0.355976},
      {"a2V_BKstar","a0T1_BKstar",-0.354268},
      {"a2V_BKstar","a1T1_BKstar",-0.0783557},
      {"a2V_BKstar","a2T1_BKstar",0.338328},
      {"a2V_BKstar","a0T2_BKstar",-0.352355},
      {"a2V_BKstar","a1T2_BKstar",-0.0541063},
      {"a2V_BKstar","a2T2_BKstar",0.38493},
      {"a2V_BKstar","a0T23_BKstar",0.016202},
      {"a2V_BKstar","a1T23_BKstar",6.71602e-05},
      {"a2V_BKstar","a2T23_BKstar",0.0569126},
      {"a0T1_BKstar","a1T1_BKstar",0.713877},
      {"a0T1_BKstar","a2T1_BKstar",-0.473872},
      {"a0T1_BKstar","a0T2_BKstar",0.999},
      {"a0T1_BKstar","a1T2_BKstar",0.686349},
      {"a0T1_BKstar","a2T2_BKstar",-0.229853},
      {"a0T1_BKstar","a0T23_BKstar",0.515437},
      {"a0T1_BKstar","a1T23_BKstar",0.236189},
      {"a0T1_BKstar","a2T23_BKstar",-0.310316},
      {"a1T1_BKstar","a2T1_BKstar",0.00093065},
      {"a1T1_BKstar","a0T2_BKstar",0.706907},
      {"a1T1_BKstar","a1T2_BKstar",0.881232},
      {"a1T1_BKstar","a2T2_BKstar",0.365115},
      {"a1T1_BKstar","a0T23_BKstar",0.601021},
      {"a1T1_BKstar","a1T23_BKstar",0.427022},
      {"a1T1_BKstar","a2T23_BKstar",-0.210129},
      {"a2T1_BKstar","a0T2_BKstar",-0.48143},
      {"a2T1_BKstar","a1T2_BKstar",-0.239308},
      {"a2T1_BKstar","a2T2_BKstar",0.626511},
      {"a2T1_BKstar","a0T23_BKstar",-0.27088},
      {"a2T1_BKstar","a1T23_BKstar",-0.105814},
      {"a2T1_BKstar","a2T23_BKstar",0.447129},
      {"a0T2_BKstar","a1T2_BKstar",0.692809},
      {"a0T2_BKstar","a2T2_BKstar",-0.227854},
      {"a0T2_BKstar","a0T23_BKstar",0.517204},
      {"a0T2_BKstar","a1T23_BKstar",0.237634},
      {"a0T2_BKstar","a2T23_BKstar",-0.313381},
      {"a1T2_BKstar","a2T2_BKstar",0.382807},
      {"a1T2_BKstar","a0T23_BKstar",0.658722},
      {"a1T2_BKstar","a1T23_BKstar",0.46875},
      {"a1T2_BKstar","a2T23_BKstar",-0.326765},
      {"a2T2_BKstar","a0T23_BKstar",0.0969361},
      {"a2T2_BKstar","a1T23_BKstar",0.216438},
      {"a2T2_BKstar","a2T23_BKstar",0.27039},
      {"a0T23_BKstar","a1T23_BKstar",0.327243},
      {"a0T23_BKstar","a2T23_BKstar",-0.711856},
      {"a1T23_BKstar","a2T23_BKstar",0.287454},
      {"a0A0_Bsphi","a1A0_Bsphi",0.687348},
      {"a0A0_Bsphi","a2A0_Bsphi",-0.349657},
      {"a0A0_Bsphi","a0A1_Bsphi",0.206736},
      {"a0A0_Bsphi","a1A1_Bsphi",-0.554941},
      {"a0A0_Bsphi","a2A1_Bsphi",-0.597622},
      {"a0A0_Bsphi","a0A12_Bsphi",0.996},
      {"a0A0_Bsphi","a1A12_Bsphi",0.772589},
      {"a0A0_Bsphi","a2A12_Bsphi",0.42463},
      {"a0A0_Bsphi","a0V_Bsphi",0.187815},
      {"a0A0_Bsphi","a1V_Bsphi",-0.422447},
      {"a0A0_Bsphi","a2V_Bsphi",-0.456553},
      {"a0A0_Bsphi","a0T1_Bsphi",0.0865627},
      {"a0A0_Bsphi","a1T1_Bsphi",-0.297008},
      {"a0A0_Bsphi","a2T1_Bsphi",0.0671143},
      {"a0A0_Bsphi","a0T2_Bsphi",0.0885951},
      {"a0A0_Bsphi","a1T2_Bsphi",-0.457242},
      {"a0A0_Bsphi","a2T2_Bsphi",-0.196984},
      {"a0A0_Bsphi","a0T23_Bsphi",0.703702},
      {"a0A0_Bsphi","a1T23_Bsphi",0.728003},
      {"a0A0_Bsphi","a2T23_Bsphi",0.602114},
      {"a1A0_Bsphi","a2A0_Bsphi",0.121069},
      {"a1A0_Bsphi","a0A1_Bsphi",0.0064726},
      {"a1A0_Bsphi","a1A1_Bsphi",-0.609677},
      {"a1A0_Bsphi","a2A1_Bsphi",-0.668983},
      {"a1A0_Bsphi","a0A12_Bsphi",0.668299},
      {"a1A0_Bsphi","a1A12_Bsphi",0.791667},
      {"a1A0_Bsphi","a2A12_Bsphi",0.646183},
      {"a1A0_Bsphi","a0V_Bsphi",0.0174825},
      {"a1A0_Bsphi","a1V_Bsphi",-0.52219},
      {"a1A0_Bsphi","a2V_Bsphi",-0.662488},
      {"a1A0_Bsphi","a0T1_Bsphi",0.0545029},
      {"a1A0_Bsphi","a1T1_Bsphi",-0.247063},
      {"a1A0_Bsphi","a2T1_Bsphi",-0.0128991},
      {"a1A0_Bsphi","a0T2_Bsphi",0.0486156},
      {"a1A0_Bsphi","a1T2_Bsphi",-0.525849},
      {"a1A0_Bsphi","a2T2_Bsphi",-0.32988},
      {"a1A0_Bsphi","a0T23_Bsphi",0.640094},
      {"a1A0_Bsphi","a1T23_Bsphi",0.775898},
      {"a1A0_Bsphi","a2T23_Bsphi",0.645661},
      {"a2A0_Bsphi","a0A1_Bsphi",0.091877},
      {"a2A0_Bsphi","a1A1_Bsphi",0.0390177},
      {"a2A0_Bsphi","a2A1_Bsphi",0.272613},
      {"a2A0_Bsphi","a0A12_Bsphi",-0.367399},
      {"a2A0_Bsphi","a1A12_Bsphi",-0.290877},
      {"a2A0_Bsphi","a2A12_Bsphi",0.166971},
      {"a2A0_Bsphi","a0V_Bsphi",0.207248},
      {"a2A0_Bsphi","a1V_Bsphi",-0.0243012},
      {"a2A0_Bsphi","a2V_Bsphi",-0.0754036},
      {"a2A0_Bsphi","a0T1_Bsphi",0.182472},
      {"a2A0_Bsphi","a1T1_Bsphi",-0.00649649},
      {"a2A0_Bsphi","a2T1_Bsphi",-0.220386},
      {"a2A0_Bsphi","a0T2_Bsphi",0.174415},
      {"a2A0_Bsphi","a1T2_Bsphi",-0.0908376},
      {"a2A0_Bsphi","a2T2_Bsphi",-0.260704},
      {"a2A0_Bsphi","a0T23_Bsphi",-0.255547},
      {"a2A0_Bsphi","a1T23_Bsphi",-0.197701},
      {"a2A0_Bsphi","a2T23_Bsphi",-0.0842799},
      {"a0A1_Bsphi","a1A1_Bsphi",0.250919},
      {"a0A1_Bsphi","a2A1_Bsphi",-0.105449},
      {"a0A1_Bsphi","a0A12_Bsphi",0.200928},
      {"a0A1_Bsphi","a1A12_Bsphi",-0.13586},
      {"a0A1_Bsphi","a2A12_Bsphi",-0.146239},
      {"a0A1_Bsphi","a0V_Bsphi",0.406354},
      {"a0A1_Bsphi","a1V_Bsphi",0.180876},
      {"a0A1_Bsphi","a2V_Bsphi",0.109717},
      {"a0A1_Bsphi","a0T1_Bsphi",0.267591},
      {"a0A1_Bsphi","a1T1_Bsphi",0.101937},
      {"a0A1_Bsphi","a2T1_Bsphi",-0.164558},
      {"a0A1_Bsphi","a0T2_Bsphi",0.275557},
      {"a0A1_Bsphi","a1T2_Bsphi",0.190196},
      {"a0A1_Bsphi","a2T2_Bsphi",-0.0562319},
      {"a0A1_Bsphi","a0T23_Bsphi",0.0931432},
      {"a0A1_Bsphi","a1T23_Bsphi",-0.0724199},
      {"a0A1_Bsphi","a2T23_Bsphi",-0.155631},
      {"a1A1_Bsphi","a2A1_Bsphi",0.78406},
      {"a1A1_Bsphi","a0A12_Bsphi",-0.549042},
      {"a1A1_Bsphi","a1A12_Bsphi",-0.662113},
      {"a1A1_Bsphi","a2A12_Bsphi",-0.589264},
      {"a1A1_Bsphi","a0V_Bsphi",-0.00746984},
      {"a1A1_Bsphi","a1V_Bsphi",0.642019},
      {"a1A1_Bsphi","a2V_Bsphi",0.552109},
      {"a1A1_Bsphi","a0T1_Bsphi",0.142301},
      {"a1A1_Bsphi","a1T1_Bsphi",0.571088},
      {"a1A1_Bsphi","a2T1_Bsphi",-0.22105},
      {"a1A1_Bsphi","a0T2_Bsphi",0.139241},
      {"a1A1_Bsphi","a1T2_Bsphi",0.794425},
      {"a1A1_Bsphi","a2T2_Bsphi",0.195856},
      {"a1A1_Bsphi","a0T23_Bsphi",-0.601549},
      {"a1A1_Bsphi","a1T23_Bsphi",-0.673956},
      {"a1A1_Bsphi","a2T23_Bsphi",-0.681441},
      {"a2A1_Bsphi","a0A12_Bsphi",-0.58909},
      {"a2A1_Bsphi","a1A12_Bsphi",-0.716176},
      {"a2A1_Bsphi","a2A12_Bsphi",-0.453332},
      {"a2A1_Bsphi","a0V_Bsphi",0.0402144},
      {"a2A1_Bsphi","a1V_Bsphi",0.47354},
      {"a2A1_Bsphi","a2V_Bsphi",0.497547},
      {"a2A1_Bsphi","a0T1_Bsphi",0.156373},
      {"a2A1_Bsphi","a1T1_Bsphi",0.388671},
      {"a2A1_Bsphi","a2T1_Bsphi",-0.169048},
      {"a2A1_Bsphi","a0T2_Bsphi",0.148962},
      {"a2A1_Bsphi","a1T2_Bsphi",0.548463},
      {"a2A1_Bsphi","a2T2_Bsphi",0.124908},
      {"a2A1_Bsphi","a0T23_Bsphi",-0.652782},
      {"a2A1_Bsphi","a1T23_Bsphi",-0.693798},
      {"a2A1_Bsphi","a2T23_Bsphi",-0.569813},
      {"a0A12_Bsphi","a1A12_Bsphi",0.775667},
      {"a0A12_Bsphi","a2A12_Bsphi",0.424398},
      {"a0A12_Bsphi","a0V_Bsphi",0.184513},
      {"a0A12_Bsphi","a1V_Bsphi",-0.420811},
      {"a0A12_Bsphi","a2V_Bsphi",-0.447059},
      {"a0A12_Bsphi","a0T1_Bsphi",0.0773811},
      {"a0A12_Bsphi","a1T1_Bsphi",-0.304431},
      {"a0A12_Bsphi","a2T1_Bsphi",0.0722582},
      {"a0A12_Bsphi","a0T2_Bsphi",0.0795224},
      {"a0A12_Bsphi","a1T2_Bsphi",-0.455068},
      {"a0A12_Bsphi","a2T2_Bsphi",-0.187699},
      {"a0A12_Bsphi","a0T23_Bsphi",0.698961},
      {"a0A12_Bsphi","a1T23_Bsphi",0.723823},
      {"a0A12_Bsphi","a2T23_Bsphi",0.599831},
      {"a1A12_Bsphi","a2A12_Bsphi",0.803423},
      {"a1A12_Bsphi","a0V_Bsphi",-0.0661659},
      {"a1A12_Bsphi","a1V_Bsphi",-0.502237},
      {"a1A12_Bsphi","a2V_Bsphi",-0.553914},
      {"a1A12_Bsphi","a0T1_Bsphi",-0.148005},
      {"a1A12_Bsphi","a1T1_Bsphi",-0.406181},
      {"a1A12_Bsphi","a2T1_Bsphi",0.0686011},
      {"a1A12_Bsphi","a0T2_Bsphi",-0.149814},
      {"a1A12_Bsphi","a1T2_Bsphi",-0.597757},
      {"a1A12_Bsphi","a2T2_Bsphi",-0.236414},
      {"a1A12_Bsphi","a0T23_Bsphi",0.629096},
      {"a1A12_Bsphi","a1T23_Bsphi",0.8456},
      {"a1A12_Bsphi","a2T23_Bsphi",0.733289},
      {"a2A12_Bsphi","a0V_Bsphi",0.0196031},
      {"a2A12_Bsphi","a1V_Bsphi",-0.450477},
      {"a2A12_Bsphi","a2V_Bsphi",-0.475719},
      {"a2A12_Bsphi","a0T1_Bsphi",-0.124325},
      {"a2A12_Bsphi","a1T1_Bsphi",-0.459704},
      {"a2A12_Bsphi","a2T1_Bsphi",-0.0415268},
      {"a2A12_Bsphi","a0T2_Bsphi",-0.129832},
      {"a2A12_Bsphi","a1T2_Bsphi",-0.651114},
      {"a2A12_Bsphi","a2T2_Bsphi",-0.338525},
      {"a2A12_Bsphi","a0T23_Bsphi",0.369588},
      {"a2A12_Bsphi","a1T23_Bsphi",0.656874},
      {"a2A12_Bsphi","a2T23_Bsphi",0.656771},
      {"a0V_Bsphi","a1V_Bsphi",0.247761},
      {"a0V_Bsphi","a2V_Bsphi",-0.123124},
      {"a0V_Bsphi","a0T1_Bsphi",0.526268},
      {"a0V_Bsphi","a1T1_Bsphi",-0.0393931},
      {"a0V_Bsphi","a2T1_Bsphi",-0.415068},
      {"a0V_Bsphi","a0T2_Bsphi",0.535066},
      {"a0V_Bsphi","a1T2_Bsphi",-0.0392003},
      {"a0V_Bsphi","a2T2_Bsphi",-0.368212},
      {"a0V_Bsphi","a0T23_Bsphi",0.160071},
      {"a0V_Bsphi","a1T23_Bsphi",0.0458389},
      {"a0V_Bsphi","a2T23_Bsphi",0.0186265},
      {"a1V_Bsphi","a2V_Bsphi",0.739804},
      {"a1V_Bsphi","a0T1_Bsphi",0.218671},
      {"a1V_Bsphi","a1T1_Bsphi",0.571275},
      {"a1V_Bsphi","a2T1_Bsphi",-0.444945},
      {"a1V_Bsphi","a0T2_Bsphi",0.223869},
      {"a1V_Bsphi","a1T2_Bsphi",0.633854},
      {"a1V_Bsphi","a2T2_Bsphi",-0.0725463},
      {"a1V_Bsphi","a0T23_Bsphi",-0.617296},
      {"a1V_Bsphi","a1T23_Bsphi",-0.617176},
      {"a1V_Bsphi","a2T23_Bsphi",-0.629361},
      {"a2V_Bsphi","a0T1_Bsphi",-0.099119},
      {"a2V_Bsphi","a1T1_Bsphi",0.323867},
      {"a2V_Bsphi","a2T1_Bsphi",-0.0726924},
      {"a2V_Bsphi","a0T2_Bsphi",-0.0931479},
      {"a2V_Bsphi","a1T2_Bsphi",0.50992},
      {"a2V_Bsphi","a2T2_Bsphi",0.213584},
      {"a2V_Bsphi","a0T23_Bsphi",-0.677613},
      {"a2V_Bsphi","a1T23_Bsphi",-0.651525},
      {"a2V_Bsphi","a2T23_Bsphi",-0.547041},
      {"a0T1_Bsphi","a1T1_Bsphi",0.448903},
      {"a0T1_Bsphi","a2T1_Bsphi",-0.56142},
      {"a0T1_Bsphi","a0T2_Bsphi",0.996},
      {"a0T1_Bsphi","a1T2_Bsphi",0.189067},
      {"a0T1_Bsphi","a2T2_Bsphi",-0.51613},
      {"a0T1_Bsphi","a0T23_Bsphi",0.106616},
      {"a0T1_Bsphi","a1T23_Bsphi",-0.176712},
      {"a0T1_Bsphi","a2T23_Bsphi",-0.210211},
      {"a1T1_Bsphi","a2T1_Bsphi",-0.169757},
      {"a1T1_Bsphi","a0T2_Bsphi",0.423782},
      {"a1T1_Bsphi","a1T2_Bsphi",0.660917},
      {"a1T1_Bsphi","a2T2_Bsphi",-0.00861968},
      {"a1T1_Bsphi","a0T23_Bsphi",-0.321698},
      {"a1T1_Bsphi","a1T23_Bsphi",-0.487256},
      {"a1T1_Bsphi","a2T23_Bsphi",-0.54514},
      {"a2T1_Bsphi","a0T2_Bsphi",-0.561865},
      {"a2T1_Bsphi","a1T2_Bsphi",-0.137407},
      {"a2T1_Bsphi","a2T2_Bsphi",0.711765},
      {"a2T1_Bsphi","a0T23_Bsphi",0.220799},
      {"a2T1_Bsphi","a1T23_Bsphi",0.232593},
      {"a2T1_Bsphi","a2T23_Bsphi",0.364997},
      {"a0T2_Bsphi","a1T2_Bsphi",0.212618},
      {"a0T2_Bsphi","a2T2_Bsphi",-0.495299},
      {"a0T2_Bsphi","a0T23_Bsphi",0.109431},
      {"a0T2_Bsphi","a1T23_Bsphi",-0.176517},
      {"a0T2_Bsphi","a2T23_Bsphi",-0.208143},
      {"a1T2_Bsphi","a2T2_Bsphi",0.410351},
      {"a1T2_Bsphi","a0T23_Bsphi",-0.455203},
      {"a1T2_Bsphi","a1T23_Bsphi",-0.62091},
      {"a1T2_Bsphi","a2T23_Bsphi",-0.661853},
      {"a2T2_Bsphi","a0T23_Bsphi",-0.0306923},
      {"a2T2_Bsphi","a1T23_Bsphi",-0.114741},
      {"a2T2_Bsphi","a2T23_Bsphi",0.00905559},
      {"a0T23_Bsphi","a1T23_Bsphi",0.777787},
      {"a0T23_Bsphi","a2T23_Bsphi",0.573289},
      {"a1T23_Bsphi","a2T23_Bsphi",0.899178},
      };

      // Reserve memory
      double **corr=(double  **) malloc((nNuisance+1)*sizeof(double *));  // Nuisance parameter correlations
      for(int iObservable = 0; iObservable <= nNuisance; ++iObservable) {
          corr[iObservable]=(double *) malloc((nNuisance+1)*sizeof(double));
      }
      // --- Needed for SuperIso backend

      BEreq::convert_correlation((nuiscorr *)corrnuis, byVal(ncorrnuis), (double **)corr, (char **)namenuisance, byVal(nNuisance));

      double **res;

      BEreq::get_th_covariance_nuisance(&res, (char**)obsnames, &nObservables, &param, &nuislist, (double **)corr);

      for(int iObservable=0; iObservable < nObservables; ++iObservable) {
          for(int jObservable = 0; jObservable < nObservables; ++jObservable) {
              result[obslist[iObservable]][obslist[jObservable]] = res[iObservable][jObservable];
          }
      }

      // Free memory
      for(int iObservable = 0; iObservable <= nNuisance; ++iObservable) {
          free(corr[iObservable]);
      }
      free(corr);

      if (flav_debug) {
          for(int iObservable=0; iObservable < nObservables; ++iObservable) {
              for(int jObservable = iObservable; jObservable < nObservables; ++jObservable) {
                  printf("Covariance %s - %s: %.4e\n",
                          obsnames[iObservable], obsnames[jObservable], result[obslist[iObservable]][obslist[jObservable]]);
              }
          }
      }

      if (flav_debug) {
          cout<<"Finished SI_theory_covariance"<<endl;
      }

    }

    /// Br b-> s gamma decays
    void SI_bsgamma(double &result)
    {
      using namespace Pipes::SI_bsgamma;
      if (flav_debug) cout<<"Starting SI_bsgamma"<<endl;

      parameters const& param = *Dep::SuperIso_modelinfo;
      double E_cut=1.6;
      result=BEreq::bsgamma_CONV(&param, byVal(E_cut));

      if (flav_debug) printf("BR(b->s gamma)=%.3e\n",result);
      if (flav_debug) cout<<"Finished SI_bsgamma"<<endl;
    }


    /// Br Bs->mumu decays for the untagged case (CP-averaged)
    void SI_Bsmumu_untag(double &result)
    {
      using namespace Pipes::SI_Bsmumu_untag;
      if (flav_debug) cout<<"Starting SI_Bsmumu_untag"<<endl;

      parameters const& param = *Dep::SuperIso_modelinfo;
      int flav=2;
      result=BEreq::Bsll_untag_CONV(&param, byVal(flav));

      if (flav_debug) printf("BR(Bs->mumu)_untag=%.3e\n",result);
      if (flav_debug) cout<<"Finished SI_Bsmumu_untag"<<endl;
    }


    /// Br Bs->ee decays for the untagged case (CP-averaged)
    void SI_Bsee_untag(double &result)
    {
      using namespace Pipes::SI_Bsee_untag;
      if (flav_debug) cout<<"Starting SI_Bsee_untag"<<endl;

      parameters const& param = *Dep::SuperIso_modelinfo;
      int flav=1;
      result=BEreq::Bsll_untag_CONV(&param, byVal(flav));

      if (flav_debug) printf("BR(Bs->ee)_untag=%.3e\n",result);
      if (flav_debug) cout<<"Finished SI_Bsee_untag"<<endl;
    }


    /// Br B0->mumu decays
    void SI_Bmumu(double &result)
    {
      using namespace Pipes::SI_Bmumu;
      if (flav_debug) cout<<"Starting SI_Bmumu"<<endl;

      parameters const& param = *Dep::SuperIso_modelinfo;
      int flav=2;
      result=BEreq::Bll_CONV(&param, byVal(flav));

      if (flav_debug) printf("BR(B->mumu)=%.3e\n",result);
      if (flav_debug) cout<<"Finished SI_Bmumu"<<endl;
    }


    /// Br B->tau nu_tau decays
    void SI_Btaunu(double &result)
    {
      using namespace Pipes::SI_Btaunu;
      if (flav_debug) cout<<"Starting SI_Btaunu"<<endl;

      parameters const& param = *Dep::SuperIso_modelinfo;
      result = BEreq::Btaunu(&param);

      if (flav_debug) printf("BR(B->tau nu)=%.3e\n",result);
      if (flav_debug) cout<<"Finished SI_Btaunu"<<endl;
    }


    /// Br B->D_s tau nu
    void SI_Dstaunu(double &result)
    {
      using namespace Pipes::SI_Dstaunu;
      if (flav_debug) cout<<"Starting SI_Dstaunu"<<endl;

      parameters const& param = *Dep::SuperIso_modelinfo;
      result = BEreq::Dstaunu(&param);

      if (flav_debug) printf("BR(Ds->tau nu)=%.3e\n",result);
      if (flav_debug) cout<<"Finished SI_Dstaunu"<<endl;
    }


    /// Br B->D_s mu nu
    void SI_Dsmunu(double &result)
    {
      using namespace Pipes::SI_Dsmunu;
      if (flav_debug) cout<<"Starting SI_Dsmunu"<<endl;

      parameters const& param = *Dep::SuperIso_modelinfo;
      result = BEreq::Dsmunu(&param);

      if (flav_debug) printf("BR(Ds->mu nu)=%.3e\n",result);
      if (flav_debug) cout<<"Finished SI_Dsmunu"<<endl;
    }


    /// Br D -> mu nu
    void SI_Dmunu(double &result)
    {
      using namespace Pipes::SI_Dmunu;
      if (flav_debug) cout<<"Starting SI_Dmunu"<<endl;

      parameters const& param = *Dep::SuperIso_modelinfo;
      result = BEreq::Dmunu(&param);

      if (flav_debug) printf("BR(D->mu nu)=%.3e\n",result);
      if (flav_debug) cout<<"Finished SI_Dmunu"<<endl;
    }


    /// Br B -> D tau nu
    void SI_BDtaunu(double &result)
    {
      using namespace Pipes::SI_BDtaunu;
      if (flav_debug) cout<<"Starting SI_BDtaunu"<<endl;

      parameters const& param = *Dep::SuperIso_modelinfo;
      if (param.model < 0) FlavBit_error().raise(LOCAL_INFO, "Unsupported model.");

      double q2_min_tau_D  = 3.16; // 1.776**2
      double q2_max_tau_D  = 11.6;   // (5.28-1.869)**2
      int gen_tau_D        = 3;
      int charge_tau_D     = 0;// D* is the charged version
      double obs_tau_D[3];

      result=BEreq::BRBDlnu(byVal(gen_tau_D), byVal( charge_tau_D), byVal(q2_min_tau_D), byVal(q2_max_tau_D), byVal(obs_tau_D), &param);

      if (flav_debug) printf("BR(B-> D tau nu )=%.3e\n",result);
      if (flav_debug) cout<<"Finished SI_BDtaunu"<<endl;
    }


    /// Br B -> D mu nu
    void SI_BDmunu(double &result)
    {
      using namespace Pipes::SI_BDmunu;
      if (flav_debug) cout<<"Starting SI_BDmunu"<<endl;

      parameters const& param = *Dep::SuperIso_modelinfo;
      if (param.model < 0) FlavBit_error().raise(LOCAL_INFO, "Unsupported model.");

      double q2_min_mu_D=  0.012; // 0.105*0.105
      double q2_max_mu_D=  11.6;   // (5.28-1.869)**2
      int gen_mu_D        =2;
      int charge_mu_D     =0;// D* is the charged version
      double obs_mu_D[3];

      result= BEreq::BRBDlnu(byVal(gen_mu_D), byVal( charge_mu_D), byVal(q2_min_mu_D), byVal(q2_max_mu_D), byVal(obs_mu_D), &param);

      if (flav_debug) printf("BR(B->D mu nu)=%.3e\n",result);
      if (flav_debug) cout<<"Finished SI_BDmunu"<<endl;
    }


    /// Br B -> D* tau nu
    void SI_BDstartaunu(double &result)
    {
      using namespace Pipes::SI_BDstartaunu;
      if (flav_debug) cout<<"Starting SI_BDstartaunu"<<endl;

      parameters const& param = *Dep::SuperIso_modelinfo;
      if (param.model < 0) FlavBit_error().raise(LOCAL_INFO, "Unsupported model.");

      double q2_min_tau_Dstar = 3.16; // 1.776**2
      double q2_max_tau_Dstar = 10.67;   //(5.279-2.01027)*(5.279-2.01027);
      int gen_tau_Dstar        =3;
      int charge_tau_Dstar     =1;// D* is the charged version
      double obs_tau_Dstar[3];

      result= BEreq::BRBDstarlnu(byVal(gen_tau_Dstar), byVal( charge_tau_Dstar), byVal(q2_min_tau_Dstar), byVal(q2_max_tau_Dstar), byVal(obs_tau_Dstar), &param);

      if (flav_debug) printf("BR(B->Dstar tau nu)=%.3e\n",result);
      if (flav_debug) cout<<"Finished SI_BDstartaunu"<<endl;
    }


    /// Br B -> D* mu nu
    void SI_BDstarmunu(double &result)
    {
      using namespace Pipes::SI_BDstarmunu;
      if (flav_debug) cout<<"Starting SI_BDstarmunu"<<endl;

      parameters const& param = *Dep::SuperIso_modelinfo;
      if (param.model < 0) FlavBit_error().raise(LOCAL_INFO, "Unsupported model.");

      double q2_min_mu_Dstar = 0.012; // 0.105*0.105
      double q2_max_mu_Dstar = 10.67;   //(5.279-2.01027)*(5.279-2.01027);
      int gen_mu_Dstar        =2;
      int charge_mu_Dstar     =1;// D* is the charged version
      double obs_mu_Dstar[3];

      result=BEreq::BRBDstarlnu(byVal(gen_mu_Dstar), byVal( charge_mu_Dstar), byVal(q2_min_mu_Dstar), byVal(q2_max_mu_Dstar), byVal(obs_mu_Dstar), &param);

      if (flav_debug) printf("BR(B->Dstar mu nu)=%.3e\n",result);
      if (flav_debug) cout<<"Finished SI_BDstarmunu"<<endl;
    }


    ///  B-> D tau nu / B-> D e nu decays
    void SI_RD(double &result)
    {
      using namespace Pipes::SI_RD;
      if (flav_debug) cout<<"Starting SI_RD"<<endl;

      parameters const& param = *Dep::SuperIso_modelinfo;
      result = BEreq::BDtaunu_BDenu(&param);

      if (flav_debug) printf("BR(B->D tau nu)/BR(B->D e nu)=%.3e\n",result);
      if (flav_debug) cout<<"Finished SI_RD"<<endl;
    }


    ///  B->D* tau nu / B-> D* e nu decays
    void SI_RDstar(double &result)
    {
      using namespace Pipes::SI_RDstar;
      if (flav_debug) cout<<"Starting SI_RDstart"<<endl;

      parameters const& param = *Dep::SuperIso_modelinfo;
      result = BEreq::BDstartaunu_BDstarenu(&param);

      if (flav_debug) printf("BR(B->D* tau nu)/BR(B->D* e nu)=%.3e\n",result);
      if (flav_debug) cout<<"Finished SI_RD*"<<endl;
    }


    /// B->K mu nu / B-> pi mu nu
    void SI_Rmu(double &result)
    {
      using namespace Pipes::SI_Rmu;
      if (flav_debug) cout<<"Starting SI_Rmu"<<endl;

      parameters const& param = *Dep::SuperIso_modelinfo;
      result = BEreq::Kmunu_pimunu(&param);

      if (flav_debug) printf("R_mu=BR(K->mu nu)/BR(pi->mu nu)=%.3e\n",result);
      if (flav_debug) cout<<"Finished SI_Rmu"<<endl;
    }


    /// 2-to-3-body decay ratio for semileptonic K and pi decays
    void SI_Rmu23(double &result)
    {
      using namespace Pipes::SI_Rmu23;
      if (flav_debug) cout<<"Starting SI_Rmu23"<<endl;

      parameters const& param = *Dep::SuperIso_modelinfo;
      result = BEreq::Rmu23(&param);

      if (flav_debug) printf("Rmu23=%.3e\n",result);
      if (flav_debug) cout<<"Finished SI_Rmu23"<<endl;
    }


    /// Delta_0 (CP-averaged isospin asymmetry of B -> K* gamma)
    void SI_delta0(double &result)
    {
      using namespace Pipes::SI_delta0;
      if (flav_debug) cout<<"Starting SI_delta0"<<endl;

      parameters const& param = *Dep::SuperIso_modelinfo;
      result=BEreq::delta0_CONV(&param);

      if (flav_debug) printf("Delta0(B->K* gamma)=%.3e\n",result);
      if (flav_debug) cout<<"Finished SI_delta0"<<endl;
    }


    /// Inclusive branching fraction B -> X_s mu mu at low q^2
    void SI_BRBXsmumu_lowq2(double &result)
    {
      using namespace Pipes::SI_BRBXsmumu_lowq2;
      if (flav_debug) cout<<"Starting SI_BRBXsmumu_lowq2"<<endl;

      parameters const& param = *Dep::SuperIso_modelinfo;
      result=BEreq::BRBXsmumu_lowq2_CONV(&param);

      if (flav_debug) printf("BR(B->Xs mu mu)_lowq2=%.3e\n",result);
      if (flav_debug) cout<<"Finished SI_BRBXsmumu_lowq2"<<endl;
    }


    /// Inclusive branching fraction B -> X_s mu mu at high q^2
    void SI_BRBXsmumu_highq2(double &result)
    {
      using namespace Pipes::SI_BRBXsmumu_highq2;
      if (flav_debug) cout<<"Starting SI_BRBXsmumu_highq2"<<endl;

      parameters const& param = *Dep::SuperIso_modelinfo;
      result=BEreq::BRBXsmumu_highq2_CONV(&param);

      if (flav_debug) printf("BR(B->Xs mu mu)_highq2=%.3e\n",result);
      if (flav_debug) cout<<"Finished SI_BRBXsmumu_highq2"<<endl;
    }


    /// Forward-backward asymmetry of B -> X_s mu mu at low q^2
    void SI_A_BXsmumu_lowq2(double &result)
    {
      using namespace Pipes::SI_A_BXsmumu_lowq2;
      if (flav_debug) cout<<"Starting SI_A_BXsmumu_lowq2"<<endl;

      parameters const& param = *Dep::SuperIso_modelinfo;
      result=BEreq::A_BXsmumu_lowq2_CONV(&param);

      if (flav_debug) printf("AFB(B->Xs mu mu)_lowq2=%.3e\n",result);
      if (flav_debug) cout<<"Finished SI_A_BXsmumu_lowq2"<<endl;
    }


    /// Forward-backward asymmetry of B -> X_s mu mu at high q^2
    void SI_A_BXsmumu_highq2(double &result)
    {
      using namespace Pipes::SI_A_BXsmumu_highq2;
      if (flav_debug) cout<<"Starting SI_A_BXsmumu_highq2"<<endl;

      parameters const& param = *Dep::SuperIso_modelinfo;
      result=BEreq::A_BXsmumu_highq2_CONV(&param);

      if (flav_debug) printf("AFB(B->Xs mu mu)_highq2=%.3e\n",result);
      if (flav_debug) cout<<"Finished SI_A_BXsmumu_highq2"<<endl;
    }


    /// Zero crossing of the forward-backward asymmetry of B -> X_s mu mu
    void SI_A_BXsmumu_zero(double &result)
    {
      using namespace Pipes::SI_A_BXsmumu_zero;
      if (flav_debug) cout<<"Starting SI_A_BXsmumu_zero"<<endl;

      parameters const& param = *Dep::SuperIso_modelinfo;
      result=BEreq::A_BXsmumu_zero_CONV(&param);

      if (flav_debug) printf("AFB(B->Xs mu mu)_zero=%.3e\n",result);
      if (flav_debug) cout<<"Finished SI_A_BXsmumu_zero"<<endl;
    }


    /// Inclusive branching fraction B -> X_s tau tau at high q^2
    void SI_BRBXstautau_highq2(double &result)
    {
      using namespace Pipes::SI_BRBXstautau_highq2;
      if (flav_debug) cout<<"Starting SI_BRBXstautau_highq2"<<endl;

      parameters const& param = *Dep::SuperIso_modelinfo;
      result=BEreq::BRBXstautau_highq2_CONV(&param);

      if (flav_debug) printf("BR(B->Xs tau tau)_highq2=%.3e\n",result);
      if (flav_debug) cout<<"Finished SI_BRBXstautau_highq2"<<endl;
    }


    /// Forward-backward asymmetry of B -> X_s tau tau at high q^2
    void SI_A_BXstautau_highq2(double &result)
    {
      using namespace Pipes::SI_A_BXstautau_highq2;
      if (flav_debug) cout<<"Starting SI_A_BXstautau_highq2"<<endl;

      parameters const& param = *Dep::SuperIso_modelinfo;
      result=BEreq::A_BXstautau_highq2_CONV(&param);

      if (flav_debug) printf("AFB(B->Xs tau tau)_highq2=%.3e\n",result);
      if (flav_debug) cout<<"Finished SI_A_BXstautau_highq2"<<endl;
    }


    /// B-> K* mu mu observables in different q^2 bins
    /// @{
    #define DEFINE_BKSTARMUMU(Q2MIN, Q2MAX, Q2MIN_TAG, Q2MAX_TAG)                         \
    void CAT_4(SI_BKstarmumu_,Q2MIN_TAG,_,Q2MAX_TAG)(Flav_KstarMuMu_obs &result)          \
    {                                                                                       \
      using namespace Pipes::CAT_4(SI_BKstarmumu_,Q2MIN_TAG,_,Q2MAX_TAG);                 \
      if (flav_debug) cout<<"Starting " STRINGIFY(CAT_4(SI_BKstarmumu_,Q2MIN_TAG,_,Q2MAX_TAG))<<endl; \
      parameters const& param = *Dep::SuperIso_modelinfo;                                   \
      result=BEreq::BKstarmumu_CONV(&param, Q2MIN, Q2MAX);                                \
      if (flav_debug) cout<<"Finished " STRINGIFY(CAT_4(SI_BKstarmumu_,Q2MIN_TAG,_,Q2MAX_TAG))<<endl; \
    }
    DEFINE_BKSTARMUMU(0.1, 0.98, 0p1, 0p98)
    DEFINE_BKSTARMUMU(1.1, 2.5, 11, 25)
    DEFINE_BKSTARMUMU(2.5, 4.0, 25, 40)
    DEFINE_BKSTARMUMU(4.0, 6.0, 40, 60)
    DEFINE_BKSTARMUMU(6.0, 8.0, 60, 80)
    DEFINE_BKSTARMUMU(15., 17., 15, 17)
    DEFINE_BKSTARMUMU(17., 19., 17, 19)
    DEFINE_BKSTARMUMU(15., 19., 15, 19)
    /// @}
    #undef DEFINE_BKSTARMUMU

    /// RK* in low q^2
    void SI_RKstar_0045_11(double &result)
    {
      using namespace Pipes::SI_RKstar_0045_11;
      if (flav_debug) cout<<"Starting SI_RKstar_0045_11"<<endl;

      parameters const& param = *Dep::SuperIso_modelinfo;
      result=BEreq::RKstar_CONV(&param,0.045,1.1);

      if (flav_debug) printf("RK*_lowq2=%.3e\n",result);
      if (flav_debug) cout<<"Finished SI_RKstar_0045_11"<<endl;
    }

    /// RK* in intermediate q^2
    void SI_RKstar_11_60(double &result)
    {
      using namespace Pipes::SI_RKstar_11_60;
      if (flav_debug) cout<<"Starting SI_RKstar_11_60"<<endl;

      parameters const& param = *Dep::SuperIso_modelinfo;
      result=BEreq::RKstar_CONV(&param,1.1,6.0);

      if (flav_debug) printf("RK*_intermq2=%.3e\n",result);
      if (flav_debug) cout<<"Finished SI_RKstar_11_60"<<endl;
    }

    // RK* for RHN, using same approximations as RK, low q^2
    void RHN_RKstar_0045_11(double &result)
    {
      using namespace Pipes::RHN_RKstar_0045_11;
      SMInputs sminputs = *Dep::SMINPUTS;
      Eigen::Matrix3cd Theta = *Dep::SeesawI_Theta;
      std::vector<double> mN = {*Param["M_1"],*Param["M_2"],*Param["M_3"]};
      double mt = *Param["mT"];

      if (flav_debug) cout << "Starting RHN_RKstar_0045_11" << endl;

      const double mW = sminputs.mW;
      const double sinW2 = sqrt(1.0 - pow(sminputs.mW/sminputs.mZ,2));

      // NNLL calculation of SM Wilson coefficients from 1712.01593 and 0811.1214
      const double C10_SM = -4.103;
      const double C9_SM = 4.211;

      // Wilson coefficients for the RHN model, from 1706.07570
      std::complex<double> C10_mu = {0.0, 0.0}, C10_e = {0.0, 0.0};
      for(int i=0; i<3; i++)
      {
        C10_mu += 1.0/(4.0*sinW2)*Theta.adjoint()(i,1)*Theta(1,i) * LoopFunctions::E(pow(mt/mW,2),pow(mN[i]/mW,2));
        C10_e += 1.0/(4.0*sinW2)*Theta.adjoint()(i,0)*Theta(0,i) * LoopFunctions::E(pow(mt/mW,2),pow(mN[i]/mW,2));
      }
      std::complex<double> C9_mu = - C10_mu, C9_e = -C10_e;

      // Aproximated solution from eq A.3 in 1408.4097
      result =  std::norm(C10_SM + C10_mu) + std::norm(C9_SM + C9_mu);
      result /= std::norm(C10_SM + C10_e) + std::norm(C9_SM + C9_e);

      if (flav_debug) cout << "RK = " << result << endl;
      if (flav_debug) cout << "Finished RHN_RKstar_0045_11" << endl;

    }

    // RK* for RHN, using same approximations as RK, intermediate q^2
    void RHN_RKstar_11_60(double &result)
    {
      using namespace Pipes::RHN_RKstar_11_60;
      SMInputs sminputs = *Dep::SMINPUTS;
      Eigen::Matrix3cd Theta = *Dep::SeesawI_Theta;
      std::vector<double> mN = {*Param["M_1"],*Param["M_2"],*Param["M_3"]};
      double mt = *Param["mT"];

      if (flav_debug) cout << "Starting RHN_RKstar_11_60" << endl;

      const double mW = sminputs.mW;
      const double sinW2 = sqrt(1.0 - pow(sminputs.mW/sminputs.mZ,2));

      // NNLL calculation of SM Wilson coefficients from 1712.01593 and 0811.1214
      const double C10_SM = -4.103;
      const double C9_SM = 4.211;

      // Wilson coefficients for the RHN model, from 1706.07570
      std::complex<double> C10_mu = {0.0, 0.0}, C10_e = {0.0, 0.0};
      for(int i=0; i<3; i++)
      {
        C10_mu += 1.0/(4.0*sinW2)*Theta.adjoint()(i,1)*Theta(1,i) * LoopFunctions::E(pow(mt/mW,2),pow(mN[i]/mW,2));
        C10_e += 1.0/(4.0*sinW2)*Theta.adjoint()(i,0)*Theta(0,i) * LoopFunctions::E(pow(mt/mW,2),pow(mN[i]/mW,2));
      }
      std::complex<double> C9_mu = - C10_mu, C9_e = -C10_e;

      // Aproximated solution from eq A.3 in 1408.4097
      result =  std::norm(C10_SM + C10_mu) + std::norm(C9_SM + C9_mu);
      result /= std::norm(C10_SM + C10_e) + std::norm(C9_SM + C9_e);

      if (flav_debug) cout << "RK = " << result << endl;
      if (flav_debug) cout << "Finished RHN_RKstar_11_60" << endl;

    }

    /// RK between 1 and 6 GeV^2
    void SI_RK(double &result)
    {
      using namespace Pipes::SI_RK;
      if (flav_debug) cout<<"Starting SI_RK"<<endl;

      parameters const& param = *Dep::SuperIso_modelinfo;
      result=BEreq::RK_CONV(&param,1.0,6.0);

      if (flav_debug) printf("RK=%.3e\n",result);
      if (flav_debug) cout<<"Finished SI_RK"<<endl;
    }

    /// RK for RHN
    void RHN_RK(double &result)
    {
      using namespace Pipes::RHN_RK;
      SMInputs sminputs = *Dep::SMINPUTS;
      Eigen::Matrix3cd Theta = *Dep::SeesawI_Theta;
      std::vector<double> mN = {*Param["M_1"],*Param["M_2"],*Param["M_3"]};
      double mt = *Param["mT"];

      if (flav_debug) cout << "Starting RHN_RK" << endl;

      const double mW = sminputs.mW;
      const double sinW2 = sqrt(1.0 - pow(sminputs.mW/sminputs.mZ,2));

      // NNLL calculation of SM Wilson coefficients from 1712.01593 and 0811.1214
      const double C10_SM = -4.103;
      const double C9_SM = 4.211;

      // Wilson coefficients for the RHN model, from 1706.07570
      std::complex<double> C10_mu = {0.0, 0.0}, C10_e = {0.0, 0.0};
      for(int i=0; i<3; i++)
      {
        C10_mu += 1.0/(4.0*sinW2)*Theta.adjoint()(i,1)*Theta(1,i) * LoopFunctions::E(pow(mt/mW,2),pow(mN[i]/mW,2));
        C10_e += 1.0/(4.0*sinW2)*Theta.adjoint()(i,0)*Theta(0,i) * LoopFunctions::E(pow(mt/mW,2),pow(mN[i]/mW,2));
      }
      std::complex<double> C9_mu = - C10_mu, C9_e = -C10_e;

      // Aproximated solution from eq A.3 in 1408.4097
      result =  std::norm(C10_SM + C10_mu) + std::norm(C9_SM + C9_mu);
      result /= std::norm(C10_SM + C10_e) + std::norm(C9_SM + C9_e);

      if (flav_debug) cout << "RK = " << result << endl;
      if (flav_debug) cout << "Finished RHN_RK" << endl;
    }

    /// Isospin asymmetry of B-> K* mu mu
    void SI_AI_BKstarmumu(double &result)
    {
      using namespace Pipes::SI_AI_BKstarmumu;
      if (flav_debug) cout<<"Starting SI_AI_BKstarmumu"<<endl;

      parameters const& param = *Dep::SuperIso_modelinfo;
      result=BEreq::SI_AI_BKstarmumu_CONV(&param);

      if (flav_debug) printf("A_I(B->K* mu mu)_lowq2=%.3e\n",result);
      if (flav_debug) cout<<"Finished SI_AI_BKstarmumu"<<endl;
    }


    /// Zero crossing of isospin asymmetry of B-> K* mu mu
    void SI_AI_BKstarmumu_zero(double &result)
    {
      using namespace Pipes::SI_AI_BKstarmumu_zero;

      if (flav_debug) cout<<"Starting SI_AI_BKstarmumu_zero"<<endl;

      parameters const& param = *Dep::SuperIso_modelinfo;
      result=BEreq::SI_AI_BKstarmumu_zero_CONV(&param);

      if (flav_debug) printf("A_I(B->K* mu mu)_zero=%.3e\n",result);
      if (flav_debug) cout<<"Finished SI_AI_BKstarmumu_zero"<<endl;
    }


    /// Flavour observables from FeynHiggs: B_s mass asymmetry, Br B_s -> mu mu, Br B -> X_s gamma
    void FH_FlavourObs(fh_FlavourObs &result)
    {
      using namespace Pipes::FH_FlavourObs;

      if (flav_debug) cout<<"Starting FH_FlavourObs"<<endl;

      fh_real bsgMSSM;     // B -> Xs gamma in MSSM
      fh_real bsgSM;       // B -> Xs gamma in SM
      fh_real deltaMsMSSM; // delta Ms in MSSM
      fh_real deltaMsSM;   // delta Ms in SM
      fh_real bsmumuMSSM;  // Bs -> mu mu in MSSM
      fh_real bsmumuSM;    // Bs -> mu mu in SM

      int error = 1;
      BEreq::FHFlavour(error, bsgMSSM, bsgSM,
           deltaMsMSSM, deltaMsSM,
           bsmumuMSSM, bsmumuSM);

      fh_FlavourObs FlavourObs;
      FlavourObs.Bsg_MSSM = bsgMSSM;
      FlavourObs.Bsg_SM = bsgSM;
      FlavourObs.deltaMs_MSSM = deltaMsMSSM;
      FlavourObs.deltaMs_SM = deltaMsSM;
      FlavourObs.Bsmumu_MSSM = bsmumuMSSM;
      FlavourObs.Bsmumu_SM = bsmumuSM;

      result = FlavourObs;
      if (flav_debug) cout<<"Finished FH_FlavourObs"<<endl;
    }


    ///These functions extract observables from a FeynHiggs flavour result
    ///@{
    void FH_bsgamma(double &result)
    {
      result = Pipes::FH_bsgamma::Dep::FH_FlavourObs->Bsg_MSSM;
    }
    void FH_Bsmumu (double &result)
    {
      result = Pipes::FH_Bsmumu::Dep::FH_FlavourObs->Bsmumu_MSSM;
    }
    void FH_DeltaMs(double &result)
    {
      result = Pipes::FH_DeltaMs::Dep::FH_FlavourObs->deltaMs_MSSM;
    }
    ///@}


    /// Measurements for electroweak penguin decays
    void b2sll_measurements(predictions_measurements_covariances &pmc)
    {
      using namespace Pipes::b2sll_measurements;

      static bool first = true;
      static int n_experiments;

      if (flav_debug) cout<<"Starting b2sll_measurements function"<<endl;
      if (flav_debug and first) cout <<"Initialising Flav Reader in b2sll_measurements"<<endl;

      // Read and calculate things based on the observed data only the first time through, as none of it depends on the model parameters.
      if (first)
      {
        pmc.LL_name="b2sll_likelihood";

        Flav_reader fread(GAMBIT_DIR  "/FlavBit/data");
        fread.debug_mode(flav_debug);

        const vector<string> observablesn = {"FL", "AFB", "S3", "S4", "S5", "S7", "S8", "S9"};
        const vector<string> observablesq = {"1.1-2.5", "2.5-4", "4-6", "6-8", "15-17", "17-19"};
        vector<string> observables;
        for (unsigned i=0;i<observablesq.size();++i)
        {
          for (unsigned j=0;j<observablesn.size();++j)
          {
            observables.push_back(observablesn[j]+"_B0Kstar0mumu_"+observablesq[i]);
          }
        }

        for (unsigned i=0;i<observables.size();++i)
        {
          fread.read_yaml_measurement("flav_data.yaml", observables[i]);
        }

        fread.initialise_matrices();
        pmc.cov_exp = fread.get_exp_cov();
        pmc.value_exp = fread.get_exp_value();
        pmc.cov_th = Kstarmumu_theory_err().get_th_cov(observables);
        n_experiments = pmc.cov_th.size1();
        pmc.value_th.resize(n_experiments,1);
        pmc.dim=n_experiments;

        // We assert that the experiments and the observables are the same size
        assert(pmc.value_exp.size1() == observables.size());

        // Init out.
        first = false;
      }

	  printf("BKstarmumu_11_25->FL=%.3e\n",Dep::BKstarmumu_11_25->FL);

      pmc.value_th(0,0)=Dep::BKstarmumu_11_25->FL;
      pmc.value_th(1,0)=Dep::BKstarmumu_11_25->AFB;
      pmc.value_th(2,0)=Dep::BKstarmumu_11_25->S3;
      pmc.value_th(3,0)=Dep::BKstarmumu_11_25->S4;
      pmc.value_th(4,0)=Dep::BKstarmumu_11_25->S5;
      pmc.value_th(5,0)=Dep::BKstarmumu_11_25->S7;
      pmc.value_th(6,0)=Dep::BKstarmumu_11_25->S8;
      pmc.value_th(7,0)=Dep::BKstarmumu_11_25->S9;

      pmc.value_th(8,0)=Dep::BKstarmumu_25_40->FL;
      pmc.value_th(9,0)=Dep::BKstarmumu_25_40->AFB;
      pmc.value_th(10,0)=Dep::BKstarmumu_25_40->S3;
      pmc.value_th(11,0)=Dep::BKstarmumu_25_40->S4;
      pmc.value_th(12,0)=Dep::BKstarmumu_25_40->S5;
      pmc.value_th(13,0)=Dep::BKstarmumu_25_40->S7;
      pmc.value_th(14,0)=Dep::BKstarmumu_25_40->S8;
      pmc.value_th(15,0)=Dep::BKstarmumu_25_40->S9;

      pmc.value_th(16,0)=Dep::BKstarmumu_40_60->FL;
      pmc.value_th(17,0)=Dep::BKstarmumu_40_60->AFB;
      pmc.value_th(18,0)=Dep::BKstarmumu_40_60->S3;
      pmc.value_th(19,0)=Dep::BKstarmumu_40_60->S4;
      pmc.value_th(20,0)=Dep::BKstarmumu_40_60->S5;
      pmc.value_th(21,0)=Dep::BKstarmumu_40_60->S7;
      pmc.value_th(22,0)=Dep::BKstarmumu_40_60->S8;
      pmc.value_th(23,0)=Dep::BKstarmumu_40_60->S9;

      pmc.value_th(24,0)=Dep::BKstarmumu_60_80->FL;
      pmc.value_th(25,0)=Dep::BKstarmumu_60_80->AFB;
      pmc.value_th(26,0)=Dep::BKstarmumu_60_80->S3;
      pmc.value_th(27,0)=Dep::BKstarmumu_60_80->S4;
      pmc.value_th(28,0)=Dep::BKstarmumu_60_80->S5;
      pmc.value_th(29,0)=Dep::BKstarmumu_60_80->S7;
      pmc.value_th(30,0)=Dep::BKstarmumu_60_80->S8;
      pmc.value_th(31,0)=Dep::BKstarmumu_60_80->S9;

      pmc.value_th(32,0)=Dep::BKstarmumu_15_17->FL;
      pmc.value_th(33,0)=Dep::BKstarmumu_15_17->AFB;
      pmc.value_th(34,0)=Dep::BKstarmumu_15_17->S3;
      pmc.value_th(35,0)=Dep::BKstarmumu_15_17->S4;
      pmc.value_th(36,0)=Dep::BKstarmumu_15_17->S5;
      pmc.value_th(37,0)=Dep::BKstarmumu_15_17->S7;
      pmc.value_th(38,0)=Dep::BKstarmumu_15_17->S8;
      pmc.value_th(39,0)=Dep::BKstarmumu_15_17->S9;

      pmc.value_th(40,0)=Dep::BKstarmumu_17_19->FL;
      pmc.value_th(41,0)=Dep::BKstarmumu_17_19->AFB;
      pmc.value_th(42,0)=Dep::BKstarmumu_17_19->S3;
      pmc.value_th(43,0)=Dep::BKstarmumu_17_19->S4;
      pmc.value_th(44,0)=Dep::BKstarmumu_17_19->S5;
      pmc.value_th(45,0)=Dep::BKstarmumu_17_19->S7;
      pmc.value_th(46,0)=Dep::BKstarmumu_17_19->S8;
      pmc.value_th(47,0)=Dep::BKstarmumu_17_19->S9;

      pmc.diff.clear();
      for (int i=0;i<n_experiments;++i)
      {
        pmc.diff.push_back(pmc.value_exp(i,0)-pmc.value_th(i,0));
      }

      if (flav_debug) cout<<"Finished b2sll_measurements function"<<endl;
    }


    /// Likelihood for electroweak penguin decays
    void b2sll_likelihood(double &result)
    {
      using namespace Pipes::b2sll_likelihood;

      if (flav_debug) cout<<"Starting b2sll_likelihood"<<endl;

      // Get experimental measurements
      predictions_measurements_covariances pmc=*Dep::b2sll_M;

      // Get experimental covariance
      boost::numeric::ublas::matrix<double> cov=pmc.cov_exp;

      // adding theory and experimenta covariance
      cov+=pmc.cov_th;

      //calculating a diff
      vector<double> diff;
      diff=pmc.diff;
      boost::numeric::ublas::matrix<double> cov_inv(pmc.dim, pmc.dim);
      InvertMatrix(cov, cov_inv);

      double Chi2=0;

      for (int i=0; i < pmc.dim; ++i)
      {
        for (int j=0; j<pmc.dim; ++j)
        {
          Chi2+= diff[i] * cov_inv(i,j)*diff[j] ;
        }
      }

      result=-0.5*Chi2;

      if (flav_debug) cout<<"Finished b2sll_likelihood"<<endl;
      if (flav_debug_LL) cout<<"Likelihood result b2sll_likelihood : "<< result<<endl;

    }


    /// Likelihood for Delta Ms
    void deltaMB_likelihood(double &result)
    {
      using namespace Pipes::deltaMB_likelihood;
      static bool th_err_absolute, first = true;
      static double exp_meas, exp_DeltaMs_err, th_err;

      if (flav_debug) cout << "Starting Delta_Ms_likelihood"<<endl;

      if (first)
      {
        Flav_reader fread(GAMBIT_DIR  "/FlavBit/data");
        fread.debug_mode(flav_debug);
        if (flav_debug) cout<<"Initialised Flav reader in Delta_Ms_likelihood"<<endl;
        fread.read_yaml_measurement("flav_data.yaml", "DeltaMs");
        fread.initialise_matrices(); // here we have a single measurement ;) so let's be sneaky:
        exp_meas = fread.get_exp_value()(0,0);
        exp_DeltaMs_err = sqrt(fread.get_exp_cov()(0,0));
        th_err = fread.get_th_err()(0,0).first;
        th_err_absolute = fread.get_th_err()(0,0).second;
        first = false;
      }

      if (flav_debug) cout << "Experiment: " << exp_meas << " " << exp_DeltaMs_err << " " << th_err << endl;

      // Now we do the stuff that actually depends on the parameters
      double theory_prediction = *Dep::DeltaMs;
      double theory_DeltaMs_err = th_err * (th_err_absolute ? 1.0 : std::abs(theory_prediction));
      if (flav_debug) cout<<"Theory prediction: "<<theory_prediction<<" +/- "<<theory_DeltaMs_err<<endl;

      /// Option profile_systematics<bool>: Use likelihood version that has been profiled over systematic errors (default false)
      bool profile = runOptions->getValueOrDef<bool>(false, "profile_systematics");

      result = Stats::gaussian_loglikelihood(theory_prediction, exp_meas, theory_DeltaMs_err, exp_DeltaMs_err, profile);
    }


    /// Likelihood for b->s gamma
    void b2sgamma_likelihood(double &result)
    {
      using namespace Pipes::b2sgamma_likelihood;

      static bool th_err_absolute, first = true;
      static double exp_meas, exp_b2sgamma_err, th_err;

      if (flav_debug) cout << "Starting b2sgamma_measurements"<<endl;

      // Read and calculate things based on the observed data only the first time through, as none of it depends on the model parameters.
      if (first)
      {
        Flav_reader fread(GAMBIT_DIR  "/FlavBit/data");
        fread.debug_mode(flav_debug);
        if (flav_debug) cout<<"Initialised Flav reader in b2sgamma_measurements"<<endl;
        fread.read_yaml_measurement("flav_data.yaml", "BR_b2sgamma");
        fread.initialise_matrices(); // here we have a single measurement ;) so let's be sneaky:
        exp_meas = fread.get_exp_value()(0,0);
        exp_b2sgamma_err = sqrt(fread.get_exp_cov()(0,0));
        th_err = fread.get_th_err()(0,0).first;
        th_err_absolute = fread.get_th_err()(0,0).second;
        first = false;
      }

      if (flav_debug) cout << "Experiment: " << exp_meas << " " << exp_b2sgamma_err << " " << th_err << endl;

      // Now we do the stuff that actually depends on the parameters
      double theory_prediction = *Dep::bsgamma;
      double theory_b2sgamma_err = th_err * (th_err_absolute ? 1.0 : std::abs(theory_prediction));
      if (flav_debug) cout<<"Theory prediction: "<<theory_prediction<<" +/- "<<theory_b2sgamma_err<<endl;

      /// Option profile_systematics<bool>: Use likelihood version that has been profiled over systematic errors (default false)
      bool profile = runOptions->getValueOrDef<bool>(false, "profile_systematics");

      result = Stats::gaussian_loglikelihood(theory_prediction, exp_meas, theory_b2sgamma_err, exp_b2sgamma_err, profile);
    }


    /// Measurements for rare purely leptonic B decays
    void b2ll_measurements(predictions_measurements_covariances &pmc)
    {
      using namespace Pipes::b2ll_measurements;

      static bool bs2mumu_err_absolute, b2mumu_err_absolute, first = true;
      static double theory_bs2mumu_err, theory_b2mumu_err;

      if (flav_debug) cout<<"Starting b2ll_measurements"<<endl;

      // Read and calculate things based on the observed data only the first time through, as none of it depends on the model parameters.
      if (first)
      {
        pmc.LL_name="b2ll_likelihood";

        Flav_reader fread(GAMBIT_DIR  "/FlavBit/data");
        fread.debug_mode(flav_debug);

        if (flav_debug) cout<<"Initiated Flav reader in b2ll_measurements"<<endl;
        fread.read_yaml_measurement("flav_data.yaml", "BR_Bs2mumu");
        fread.read_yaml_measurement("flav_data.yaml", "BR_B02mumu");
        if (flav_debug) cout<<"Finished reading b->mumu data"<<endl;

        fread.initialise_matrices();

        theory_bs2mumu_err = fread.get_th_err()(0,0).first;
        theory_b2mumu_err = fread.get_th_err()(1,0).first;
        bs2mumu_err_absolute = fread.get_th_err()(0,0).second;
        b2mumu_err_absolute = fread.get_th_err()(1,0).second;

        pmc.value_exp=fread.get_exp_value();
        pmc.cov_exp=fread.get_exp_cov();

        pmc.value_th.resize(2,1);
        pmc.cov_th.resize(2,2);

        pmc.dim=2;

        // Init over and out.
        first = false;
      }

      // Get theory prediction
      pmc.value_th(0,0)=*Dep::Bsmumu_untag;
      pmc.value_th(1,0)=*Dep::Bmumu;

      // Compute error on theory prediction and populate the covariance matrix
      double theory_bs2mumu_error = theory_bs2mumu_err * (bs2mumu_err_absolute ? 1.0 : *Dep::Bsmumu_untag);
      double theory_b2mumu_error = theory_b2mumu_err * (b2mumu_err_absolute ? 1.0 : *Dep::Bmumu);
      pmc.cov_th(0,0)=theory_bs2mumu_error*theory_bs2mumu_error;
      pmc.cov_th(0,1)=0.;
      pmc.cov_th(1,0)=0.;
      pmc.cov_th(1,1)=theory_b2mumu_error*theory_b2mumu_error;

      // Save the differences between theory and experiment
      pmc.diff.clear();
      for (int i=0;i<2;++i)
      {
        pmc.diff.push_back(pmc.value_exp(i,0)-pmc.value_th(i,0));
      }

      if (flav_debug) cout<<"Finished b2ll_measurements"<<endl;

    }


    /// Likelihood for rare purely leptonic B decays
    void b2ll_likelihood(double &result)
    {
      using namespace Pipes::b2ll_likelihood;

      if (flav_debug) cout<<"Starting b2ll_likelihood"<<endl;

      predictions_measurements_covariances pmc = *Dep::b2ll_M;

      boost::numeric::ublas::matrix<double> cov=pmc.cov_exp;

      // adding theory and experimental covariance
      cov+=pmc.cov_th;

      //calculating a diff
      vector<double> diff;
      diff=pmc.diff;

      boost::numeric::ublas::matrix<double> cov_inv(pmc.dim, pmc.dim);
      InvertMatrix(cov, cov_inv);

      // calculating the chi2
      double Chi2=0;

      for (int i=0; i < pmc.dim; ++i)
      {
        for (int j=0; j<pmc.dim; ++j)
        {
          Chi2+= diff[i] * cov_inv(i,j)*diff[j];
        }
      }

      result=-0.5*Chi2;

      if (flav_debug) cout<<"Finished b2ll_likelihood"<<endl;
      if (flav_debug_LL) cout<<"Likelihood result b2ll_likelihood : "<< result<<endl;

    }


    /// Measurements for tree-level leptonic and semileptonic B decays
    void SL_measurements(predictions_measurements_covariances &pmc)
    {
      using namespace Pipes::SL_measurements;

      const int n_experiments=8;
      static bool th_err_absolute[n_experiments], first = true;
      static double th_err[n_experiments];

      if (flav_debug) cout<<"Starting SL_measurements"<<endl;

      // Read and calculate things based on the observed data only the first time through, as none of it depends on the model parameters.
      if (first)
      {
        pmc.LL_name="SL_likelihood";

        // Read in experimental measuremens
        Flav_reader fread(GAMBIT_DIR  "/FlavBit/data");
        fread.debug_mode(flav_debug);
        if (flav_debug) cout<<"Initialised Flav reader in SL_measurements"<<endl;

        // B-> tau nu
        fread.read_yaml_measurement("flav_data.yaml", "BR_Btaunu");
        // B-> D mu nu
        fread.read_yaml_measurement("flav_data.yaml", "BR_BDmunu");
        // B-> D* mu nu
        fread.read_yaml_measurement("flav_data.yaml", "BR_BDstarmunu");
        // RD
        fread.read_yaml_measurement("flav_data.yaml", "RD");
        // RDstar
        fread.read_yaml_measurement("flav_data.yaml", "RDstar");
        // Ds-> tau nu
        fread.read_yaml_measurement("flav_data.yaml", "BR_Dstaunu");
        // Ds -> mu nu
        fread.read_yaml_measurement("flav_data.yaml", "BR_Dsmunu");
        // D -> mu nu
        fread.read_yaml_measurement("flav_data.yaml", "BR_Dmunu");

        fread.initialise_matrices();
        pmc.cov_exp=fread.get_exp_cov();
        pmc.value_exp=fread.get_exp_value();

        pmc.value_th.resize(n_experiments,1);
        // Set all entries in the covariance matrix explicitly to zero, as we will only write the diagonal ones later.
        pmc.cov_th = boost::numeric::ublas::zero_matrix<double>(n_experiments,n_experiments);
        for (int i = 0; i < n_experiments; ++i)
        {
          th_err[i] = fread.get_th_err()(i,0).first;
          th_err_absolute[i] = fread.get_th_err()(i,0).second;
        }

        pmc.dim=n_experiments;

        // Init over.
        first = false;
      }

      // R(D) is calculated assuming isospin symmetry
      double theory[8];
      // B-> tau nu SI
      theory[0] = *Dep::Btaunu;
      // B-> D mu nu
      theory[1] = *Dep::BDmunu;
      // B-> D* mu nu
      theory[2] = *Dep::BDstarmunu;
      // RD
      theory[3] = *Dep::RD;
      // RDstar
      theory[4] = *Dep::RDstar;
      // Ds-> tau nu
      theory[5] = *Dep::Dstaunu;
      // Ds -> mu nu
      theory[6] = *Dep::Dsmunu;
      // D -> mu nu
      theory[7] =*Dep::Dmunu;

      for (int i = 0; i < n_experiments; ++i)
      {
        pmc.value_th(i,0) = theory[i];
        pmc.cov_th(i,i) = th_err[i]*th_err[i] * (th_err_absolute[i] ? 1.0 : theory[i]*theory[i]);
      }
      // Add in the correlations between B-> D mu nu and RD
      pmc.cov_th(1,3) = pmc.cov_th(3,1) = -0.55 * th_err[1]*th_err[3] * (th_err_absolute[1] ? 1.0 : theory[1]) * (th_err_absolute[3] ? 1.0 : theory[3]);
      // Add in the correlations between B-> D* mu nu and RD*
      pmc.cov_th(2,4) = pmc.cov_th(4,2) = -0.62 * th_err[2]*th_err[4] * (th_err_absolute[2] ? 1.0 : theory[2]) * (th_err_absolute[4] ? 1.0 : theory[4]);

      pmc.diff.clear();
      for (int i=0;i<n_experiments;++i)
      {
        pmc.diff.push_back(pmc.value_exp(i,0)-pmc.value_th(i,0));
      }

      if (flav_debug) cout<<"Finished SL_measurements"<<endl;

    }


    /// Likelihood for tree-level leptonic and semileptonic B decays
    void SL_likelihood(double &result)
    {
      using namespace Pipes::SL_likelihood;

      if (flav_debug) cout<<"Starting SL_likelihood"<<endl;

      predictions_measurements_covariances pmc = *Dep::SL_M;

      boost::numeric::ublas::matrix<double> cov=pmc.cov_exp;

      // adding theory and experimental covariance
      cov+=pmc.cov_th;

      //calculating a diff
      vector<double> diff;
      diff=pmc.diff;

      boost::numeric::ublas::matrix<double> cov_inv(pmc.dim, pmc.dim);
      InvertMatrix(cov, cov_inv);

      double Chi2=0;
      for (int i=0; i < pmc.dim; ++i)
      {
        for (int j=0; j<pmc.dim; ++j)
        {
          Chi2+= diff[i] * cov_inv(i,j)*diff[j];
        }
      }

      result=-0.5*Chi2;

      if (flav_debug) cout<<"Finished SL_likelihood"<<endl;

      if (flav_debug_LL) cout<<"Likelihood result SL_likelihood  : "<< result<<endl;

    }

    // Helper function
    double G(const double x)
    {
      if(x)
        return (10.0 - 43.0*x + 78.0*pow(x,2) - 49.0*pow(x,3) + 4.0*pow(x,4) + 18.0*pow(x,3)*log(x)) / (3.0*pow(x - 1,4));
      else
        return 10.0/3;
    }

    // Contribution to mu -> e gamma from RHNs
    void RHN_muegamma(double &result)
    {
      using namespace Pipes::RHN_muegamma;
      SMInputs sminputs = *Dep::SMINPUTS;

      Eigen::Matrix3cd m_nu = *Dep::m_nu;
      vector<double> ml = {sminputs.mE, sminputs.mMu, sminputs.mTau};
      vector<double> mnu = {real(m_nu(0,0)), real(m_nu(1,1)), real(m_nu(2,2)), *Param["M_1"], *Param["M_2"], *Param["M_3"]};

      Eigen::Matrix3cd Theta = *Dep::SeesawI_Theta;
      Eigen::Matrix3cd Vnu = *Dep::SeesawI_Vnu;
      Eigen::Matrix<complex<double>,3,6> U;

      for(int i=0; i<3; i++)
        for(int j=0; j<3; j++)
        {
          U(i,j) = Vnu(i,j);
          U(i,j+3) = Theta(i,j);
        }

      result = pow(sminputs.mMu,5)/(4 * sminputs.alphainv);

      // Form factors
      int e = 0, mu = 1;
      complex<double> k2l = FormFactors::K2L(mu, e, sminputs, U, ml, mnu);
      complex<double> k2r = FormFactors::K2R(mu, e, sminputs, U, ml, mnu);

      result *= (norm(k2l) + norm(k2r));

      result /= Dep::mu_minus_decay_rates->width_in_GeV;

    }

    // Contribution to tau -> e gamma from RHNs
    void RHN_tauegamma(double &result)
    {
      using namespace Pipes::RHN_tauegamma;
      SMInputs sminputs = *Dep::SMINPUTS;

      Eigen::Matrix3cd m_nu = *Dep::m_nu;
      vector<double> ml = {sminputs.mE, sminputs.mMu, sminputs.mTau};
      vector<double> mnu = {real(m_nu(0,0)), real(m_nu(1,1)), real(m_nu(2,2)), *Param["M_1"], *Param["M_2"], *Param["M_3"]};

      Eigen::Matrix3cd Theta = *Dep::SeesawI_Theta;
      Eigen::Matrix3cd Vnu = *Dep::SeesawI_Vnu;
      Eigen::Matrix<complex<double>,3,6> U;

      for(int i=0; i<3; i++)
        for(int j=0; j<3; j++)
        {
          U(i,j) = Vnu(i,j);
          U(i,j+3) = Theta(i,j);
        }

      result = pow(sminputs.mTau,5)/(4*sminputs.alphainv);

      // Form factors
      int e = 0, tau = 2;
      complex<double> k2l = FormFactors::K2L(tau, e, sminputs, U, ml, mnu);
      complex<double> k2r = FormFactors::K2R(tau, e, sminputs, U, ml, mnu);

      result *= (norm(k2l) + norm(k2r));

      result /= Dep::tau_minus_decay_rates->width_in_GeV;

    }

    // Contribution to tau -> mu gamma from RHNs
    void RHN_taumugamma(double &result)
    {
      using namespace Pipes::RHN_taumugamma;
      SMInputs sminputs = *Dep::SMINPUTS;

      Eigen::Matrix3cd m_nu = *Dep::m_nu;
      vector<double> ml = {sminputs.mE, sminputs.mMu, sminputs.mTau};
      vector<double> mnu = {real(m_nu(0,0)), real(m_nu(1,1)), real(m_nu(2,2)), *Param["M_1"], *Param["M_2"], *Param["M_3"]};

      Eigen::Matrix3cd Theta = *Dep::SeesawI_Theta;
      Eigen::Matrix3cd Vnu = *Dep::SeesawI_Vnu;
      Eigen::Matrix<complex<double>,3,6> U;

      for(int i=0; i<3; i++)
        for(int j=0; j<3; j++)
        {
          U(i,j) = Vnu(i,j);
          U(i,j+3) = Theta(i,j);
        }

      result = pow(sminputs.mTau,5)/(4 * sminputs.alphainv);

      // Form factors
      int mu = 1, tau = 2;
      complex<double> k2l = FormFactors::K2L(tau, mu, sminputs, U, ml, mnu);
      complex<double> k2r = FormFactors::K2R(tau, mu, sminputs, U, ml, mnu);

      result *= (norm(k2l) + norm(k2r));

      result /= Dep::tau_minus_decay_rates->width_in_GeV;
    }

    // General contribution to l_\alpha^- -> l_\beta^- l_\gamma^- l_\delta^+ from RHNs
    double RHN_l2lll(int alpha, int beta, int gamma, int delta, SMInputs sminputs, Eigen::Matrix3cd Vnu, Eigen::Matrix3cd Theta, Eigen::Matrix3cd m_nu, double M1, double M2, double M3, double mH)
    {
      vector<double> ml = {sminputs.mE, sminputs.mMu, sminputs.mTau};
      vector<double> mnu = {real(m_nu(0,0)), real(m_nu(1,1)), real(m_nu(2,2)), M1, M2, M3};

      Eigen::Matrix<complex<double>,3,6> U;

      for(int i=0; i<3; i++)
        for(int j=0; j<3; j++)
        {
          U(i,j) = Vnu(i,j);
          U(i,j+3) = Theta(i,j);
        }

      // Form factors
      complex<double> k2l = FormFactors::K2L(alpha, beta, sminputs, U, ml, mnu);
      complex<double> k2r = FormFactors::K2R(alpha, beta, sminputs, U, ml, mnu);
      complex<double> k1r = FormFactors::K1R(alpha, beta, sminputs, U, mnu);
      complex<double> asll = FormFactors::ASLL(alpha, beta, gamma, delta, sminputs, U, ml, mnu, mH);
      complex<double> aslr = FormFactors::ASLR(alpha, beta, gamma, delta, sminputs, U, ml, mnu, mH);
      complex<double> asrl = FormFactors::ASRL(alpha, beta, gamma, delta, sminputs, U, ml, mnu, mH);
      complex<double> asrr = FormFactors::ASRR(alpha, beta, gamma, delta, sminputs, U, ml, mnu, mH);
      complex<double> avll = FormFactors::AVLL(alpha, beta, gamma, delta, sminputs, U, ml, mnu);
      complex<double> avlr = FormFactors::AVLR(alpha, beta, gamma, delta, sminputs, U, ml, mnu);
      complex<double> avrl = FormFactors::AVLL(alpha, beta, gamma, delta, sminputs, U, ml, mnu);
      complex<double> avrr = FormFactors::AVRR(alpha, beta, gamma, delta, sminputs, U, ml, mnu);

      complex<double> avhatll = avll;
      complex<double> avhatlr = avlr;
      complex<double> avhatrl = avrl + 4. * pi / sminputs.alphainv * k1r;
      complex<double> avhatrr = avrr + 4. * pi / sminputs.alphainv * k1r;

      double l2lll = 0;
      if(beta == gamma and gamma == delta) // l(alpha)- -> l(beta)- l(beta)- l(beta)+
      {
        l2lll = real(16. * pow(pi,2) / pow(sminputs.alphainv,2) * (norm(k2l) + norm(k2r)) * (16./3.*log(ml[alpha]/ml[beta]) - 22./3.) + 1./24. * (norm(asll) + norm(asrr) + 2.*norm(aslr) + 2.*norm(asrl)) + 1./3. * (2.*norm(avhatll) + 2.*norm(avhatrr) + norm(avhatlr) + norm(avhatrl)) + 4.*pi/(3.*sminputs.alphainv)*(k2l*conj(asrl - 2.*avhatrl - 4.*avhatrr) + conj(k2l)*(asrl - 2.*avhatrl - 4.*avhatrr) + k2r*conj(aslr - 2.*avhatlr - 4.*avhatll) + conj(k2r)*(aslr - 2.*avhatlr - 4.*avhatll)) - 1./6. * (aslr*conj(avhatlr) + asrl*conj(avhatrl) + conj(aslr)*avhatlr + conj(asrl)*avhatrl));
      }
      else if(gamma == delta) // l(alpha)- -> l(beta)- l(gamma)- l(gamma)+
      {
        l2lll = real(16. *pow(pi,2) / pow(sminputs.alphainv,2) * (norm(k2l) + norm(k2r)) * (16./3.*log(ml[alpha]/ml[gamma]) - 8.) + 1./12. *(norm(asll) + norm(asrr) + norm(aslr) + norm(asrl)) + 1./3. * (norm(avhatll) + norm(avhatrr) + norm(avhatlr) + norm(avhatrl)) + 8.*pi/(3.*sminputs.alphainv) * (k2l*conj(avhatrl + avhatrr) + k2r*conj(avhatlr + avhatll) + conj(k2l)*(avhatrl + avhatrr) + conj(k2r)*(avhatlr + avhatll)));
      }
      else if(beta == gamma) // l(alpha)- -> l(beta)- l(beta)- l(delta)+
      {
        l2lll = real(1./24. * (norm(asll) + norm(asrr) + 2.*norm(aslr) + 2.*norm(asrl)) + 1./3.*(2.*norm(avhatll) + 2.*norm(avhatrr) + norm(avhatlr) + norm(avhatrl)) - 1./6.*(aslr*conj(avhatlr) + asrl*conj(avhatrl) + conj(aslr)*avhatlr + conj(asrl)*avhatrl));
      }
      return l2lll;

    }

    // Contribution to mu -> e e e from RHNs
    void RHN_mueee(double &result)
    {
      using namespace Pipes::RHN_mueee;
      SMInputs sminputs = *Dep::SMINPUTS;

      Eigen::Matrix3cd m_nu = *Dep::m_nu;
      Eigen::Matrix3cd Theta = *Dep::SeesawI_Theta;
4      Eigen::Matrix3cd Vnu = *Dep::SeesawI_Vnu;

      result = pow(sminputs.mMu,5)/(512*pow(pi,3));

      int e = 0, mu = 1;
      result *=  RHN_l2lll(mu, e, e, e, sminputs, Vnu, Theta, m_nu, *Param["M_1"], *Param["M_2"], *Param["M_3"], *Param["mH"]);

      result /= Dep::mu_minus_decay_rates->width_in_GeV;

    }

    // Contribution to tau -> e e e from RHNs
    void RHN_taueee(double &result)
    {
      using namespace Pipes::RHN_taueee;
      SMInputs sminputs = *Dep::SMINPUTS;

      Eigen::Matrix3cd m_nu = *Dep::m_nu;
      Eigen::Matrix3cd Theta = *Dep::SeesawI_Theta;
      Eigen::Matrix3cd Vnu = *Dep::SeesawI_Vnu;

      result = pow(sminputs.mTau,5)/(512*pow(pi,3));

      int e = 0, tau = 2;
      result *=  RHN_l2lll(tau, e, e, e, sminputs, Vnu, Theta, m_nu, *Param["M_1"], *Param["M_2"], *Param["M_3"], *Param["mH"]);

      result /= Dep::tau_minus_decay_rates->width_in_GeV;

    }

    // Contribution to tau -> mu mu mu from RHNs
    void RHN_taumumumu(double &result)
    {
      using namespace Pipes::RHN_taumumumu;
      SMInputs sminputs = *Dep::SMINPUTS;

      Eigen::Matrix3cd m_nu = *Dep::m_nu;
      Eigen::Matrix3cd Theta = *Dep::SeesawI_Theta;
      Eigen::Matrix3cd Vnu = *Dep::SeesawI_Vnu;

      result = pow(sminputs.mTau,5)/(512*pow(pi,3));

      int mu = 1, tau = 2;
      result *=  RHN_l2lll(tau, mu, mu, mu, sminputs, Vnu, Theta, m_nu, *Param["M_1"], *Param["M_2"], *Param["M_3"], *Param["mH"]);

      result /= Dep::tau_minus_decay_rates->width_in_GeV;

    }

    // Contribution to tau^- -> mu^- e^- e^+ from RHNs
    void RHN_taumuee(double &result)
    {
      using namespace Pipes::RHN_taumuee;
      SMInputs sminputs = *Dep::SMINPUTS;

      Eigen::Matrix3cd m_nu = *Dep::m_nu;
      Eigen::Matrix3cd Theta = *Dep::SeesawI_Theta;
      Eigen::Matrix3cd Vnu = *Dep::SeesawI_Vnu;

      result = pow(sminputs.mTau,5)/(512*pow(pi,3));

      int e = 0, mu = 1, tau = 2;
      result *=  RHN_l2lll(tau, mu, e, e, sminputs, Vnu, Theta, m_nu, *Param["M_1"], *Param["M_2"], *Param["M_3"], *Param["mH"]);

      result /= Dep::tau_minus_decay_rates->width_in_GeV;
    }

    // Contribution to tau^- -> e^- e^- mu^+ from RHNs
    void RHN_taueemu(double &result)
    {
      using namespace Pipes::RHN_taueemu;
      SMInputs sminputs = *Dep::SMINPUTS;

      Eigen::Matrix3cd m_nu = *Dep::m_nu;
      Eigen::Matrix3cd Theta = *Dep::SeesawI_Theta;
      Eigen::Matrix3cd Vnu = *Dep::SeesawI_Vnu;

      result = pow(sminputs.mTau,5)/(512*pow(pi,3));

      int e = 0, mu = 1, tau = 2;
      result *=  RHN_l2lll(tau, e, e, mu, sminputs, Vnu, Theta, m_nu, *Param["M_1"], *Param["M_2"], *Param["M_3"], *Param["mH"]);

      result /= Dep::tau_minus_decay_rates->width_in_GeV;
    }

    // Contribution to tau^- -> e^- mu^- mu^+ from RHNs
    void RHN_tauemumu(double &result)
    {
      using namespace Pipes::RHN_tauemumu;
      SMInputs sminputs = *Dep::SMINPUTS;

      Eigen::Matrix3cd m_nu = *Dep::m_nu;
      Eigen::Matrix3cd Theta = *Dep::SeesawI_Theta;
      Eigen::Matrix3cd Vnu = *Dep::SeesawI_Vnu;

      result = pow(sminputs.mTau,5)/(512*pow(pi,3));

      int e = 0, mu = 1, tau = 2;
      result *=  RHN_l2lll(tau, e, mu, mu, sminputs, Vnu, Theta, m_nu, *Param["M_1"], *Param["M_2"], *Param["M_3"], *Param["mH"]);

      result /= Dep::tau_minus_decay_rates->width_in_GeV;
    }

    // Contribution to tau^- -> mu^- mu^- e^+ from RHNs
    void RHN_taumumue(double &result)
    {
      using namespace Pipes::RHN_taumumue;
      SMInputs sminputs = *Dep::SMINPUTS;

      Eigen::Matrix3cd m_nu = *Dep::m_nu;
      Eigen::Matrix3cd Theta = *Dep::SeesawI_Theta;
      Eigen::Matrix3cd Vnu = *Dep::SeesawI_Vnu;

      result = pow(sminputs.mTau,5)/(512*pow(pi,3));

      int e = 0, mu = 1, tau = 2;
      result *=  RHN_l2lll(tau, mu, mu, e, sminputs, Vnu, Theta, m_nu, *Param["M_1"], *Param["M_2"], *Param["M_3"], *Param["mH"]);

      result /= Dep::tau_minus_decay_rates->width_in_GeV;
    }

    // Contribution to mu - e conversion in Ti nucleii from RHNs
    void RHN_mueTi(double &result)
    {
      using namespace Pipes::RHN_mueTi;
      const SMInputs sminputs = *Dep::SMINPUTS;
      Eigen::Matrix3cd m_nu = *Dep::m_nu;
      Eigen::Matrix3cd Vnu = *Dep::SeesawI_Vnu;
      Eigen::Matrix3cd Theta = *Dep::SeesawI_Theta;

      vector<double> ml = {sminputs.mE, sminputs.mMu, sminputs.mTau};
      vector<double> mnu = {real(m_nu(0,0)), real(m_nu(1,1)), real(m_nu(2,2)), *Param["M_1"], *Param["M_2"], *Param["M_3"]};
      Eigen::Matrix<complex<double>,3,6> U;

      for(int i=0; i<3; i++)
        for(int j=0; j<3; j++)
        {
          U(i,j) = Vnu(i,j);
          U(i,j+3) = Theta(i,j);
        }

      // Form factors
      int e = 0, mu = 1;
      complex<double> k1r = FormFactors::K1R(mu, e, sminputs, U, mnu);
      complex<double> k2l = FormFactors::K2L(mu, e, sminputs, U, ml, mnu);
      complex<double> k2r = FormFactors::K2R(mu, e, sminputs, U, ml, mnu);

      int u = 0, d =0, s = 1;
      complex<double> CVLLu = FormFactors::CVLL(mu, e, u, u, sminputs, U, ml, mnu);
      complex<double> CVLLd = FormFactors::BVLL(mu, e, d, d, sminputs, U, ml, mnu);
      complex<double> CVLLs = FormFactors::BVLL(mu, e, s, s, sminputs, U, ml, mnu);
      complex<double> CVLRu = FormFactors::CVLR(mu, e, u, u, sminputs, U, ml, mnu);
      complex<double> CVLRd = FormFactors::BVLR(mu, e, d, d, sminputs, U, ml, mnu);
      complex<double> CVLRs = FormFactors::BVLR(mu, e, s, s, sminputs, U, ml, mnu);
      complex<double> CVRLu = FormFactors::CVRL(mu, e, u, u, sminputs, U, ml, mnu);
      complex<double> CVRLd = FormFactors::BVRL(mu, e, d, d, sminputs, U, ml, mnu);
      complex<double> CVRLs = FormFactors::BVRL(mu, e, s, s, sminputs, U, ml, mnu);
      complex<double> CVRRu = FormFactors::CVRR(mu, e, u, u, sminputs, U, ml, mnu);
      complex<double> CVRRd = FormFactors::BVRR(mu, e, d, d, sminputs, U, ml, mnu);
      complex<double> CVRRs = FormFactors::BVRR(mu, e, s, s, sminputs, U, ml, mnu);

      complex<double> CSLLu = FormFactors::CSLL(mu, e, u, u, sminputs, U, ml, mnu, *Param["mH"]);
      complex<double> CSLLd = FormFactors::BSLL(mu, e, d, d, sminputs, U, ml, mnu, *Param["mH"]);
      complex<double> CSLLs = FormFactors::BSLL(mu, e, s, s, sminputs, U, ml, mnu, *Param["mH"]);
      complex<double> CSLRu = FormFactors::CSLL(mu, e, u, u, sminputs, U, ml, mnu, *Param["mH"]);
      complex<double> CSLRd = FormFactors::BSLL(mu, e, d, d, sminputs, U, ml, mnu, *Param["mH"]);
      complex<double> CSLRs = FormFactors::BSLL(mu, e, s, s, sminputs, U, ml, mnu, *Param["mH"]);
      complex<double> CSRLu = FormFactors::CSLL(mu, e, u, u, sminputs, U, ml ,mnu, *Param["mH"]);
      complex<double> CSRLd = FormFactors::BSLL(mu, e, d, d, sminputs, U, ml ,mnu, *Param["mH"]);
      complex<double> CSRLs = FormFactors::BSLL(mu, e, s, s, sminputs, U, ml ,mnu, *Param["mH"]);
      complex<double> CSRRu = FormFactors::CSLL(mu, e, u, u, sminputs, U, ml ,mnu, *Param["mH"]);
      complex<double> CSRRd = FormFactors::BSLL(mu, e, d, d, sminputs, U, ml, mnu, *Param["mH"]);
      complex<double> CSRRs = FormFactors::BSLL(mu, e, s, s, sminputs, U, ml ,mnu, *Param["mH"]);

      double Qu = 2./3.;
      complex<double> gVLu = sqrt(2)/sminputs.GF * (4.*pi / sminputs.alphainv * Qu * (0. - k2r) - 0.5*(CVLLu + CVLRu));
      complex<double> gSLu = -1./(sqrt(2)*sminputs.GF)*(CSLLu + CSLRu);
      complex<double> gVRu = sqrt(2)/sminputs.GF * (4.*pi / sminputs.alphainv * Qu * (k1r - k2l) - 0.5*(CVRRu + CVRLu));
      complex<double> gSRu = -1./(sqrt(2)*sminputs.GF)*(CSRRu + CSRLu);

      double Qd = -1./3.;
      complex<double> gVLd = sqrt(2)/sminputs.GF * (4.*pi / sminputs.alphainv * Qd * (0. - k2r) - 0.5*(CVLLd + CVLRd));
      complex<double> gSLd = -1./(sqrt(2)*sminputs.GF)*(CSLLd + CSLRd);
      complex<double> gVRd = sqrt(2)/sminputs.GF * (4.*pi / sminputs.alphainv * Qd * (k1r - k2l) - 0.5*(CVRRd + CVRLd));
      complex<double> gSRd = -1./(sqrt(2)*sminputs.GF)*(CSRRd + CSRLd);

      double Qs = -1./3.;
      complex<double> gVLs = sqrt(2)/sminputs.GF * (4.*pi / sminputs.alphainv * Qs * (0. - k2r) - 0.5*(CVLLs + CVLRs));
      complex<double> gSLs = -1./(sqrt(2)*sminputs.GF)*(CSLLs + CSLRs);
      complex<double> gVRs = sqrt(2)/sminputs.GF * (4.*pi / sminputs.alphainv * Qs * (k1r - k2l) - 0.5*(CVRRs + CVRLs));
      complex<double> gSRs = -1./(sqrt(2)*sminputs.GF)*(CSRRs + CSRLs);

      double GVup = 2, GVdn = 2, GVdp = 1, GVun = 1, GVsp = 0, GVsn = 0;
      double GSup = 5.1, GSdn = 5.1, GSdp = 4.3, GSun = 4.3, GSsp = 2.5, GSsn = 2.5;
      complex<double> g0SL = 0.5*(gSLu*(GSup + GSun) + gSLd*(GSdp + GSdn) + gSLs*(GSsp + GSsn));
      complex<double> g0SR = 0.5*(gSRu*(GSup + GSun) + gSRd*(GSdp + GSdn) + gSRs*(GSsp + GSsn));
      complex<double> g0VL = 0.5*(gVLu*(GVup + GVun) + gVLd*(GVdp + GVdn) + gVLs*(GVsp + GVsn));
      complex<double> g0VR = 0.5*(gVRu*(GVup + GVun) + gVRd*(GVdp + GVdn) + gVRs*(GVsp + GVsn));
      complex<double> g1SL = 0.5*(gSLu*(GSup - GSun) + gSLd*(GSdp - GSdn) + gSLs*(GSsp - GSsn));
      complex<double> g1SR = 0.5*(gSRu*(GSup - GSun) + gSRd*(GSdp - GSdn) + gSRs*(GSsp - GSsn));
      complex<double> g1VL = 0.5*(gVLu*(GVup - GVun) + gVLd*(GVdp - GVdn) + gVLs*(GVsp - GVsn));
      complex<double> g1VR = 0.5*(gVRu*(GVup - GVun) + gVRd*(GVdp - GVdn) + gVRs*(GVsp - GVsn));


      // Parameters for Ti, from Table 1 in 1209.2679 for Ti
      double Z = 22, N = 26;
      double Zeff = 17.6, Fp = 0.54;
      double hbar = 6.582119514e-25; // GeV * s
      double GammaCapt = 2.59e6 * hbar;

      result = (pow(sminputs.GF,2)*pow(sminputs.mMu,5)*pow(Zeff,4)*pow(Fp,2)) / (8.*pow(pi,4)*pow(sminputs.alphainv,3)*Z*GammaCapt) * (norm((Z+N)*(g0VL + g0SL) + (Z-N)*(g1VL + g1SL)) + norm((Z+N)*(g0VR + g0SR) + (Z-N)*(g1VR + g1SR)));

    }

    // Contribution to mu - e conversion in Pb nucleii from RHNs
    void RHN_muePb(double &result)
    {
      using namespace Pipes::RHN_muePb;
      const SMInputs sminputs = *Dep::SMINPUTS;
      Eigen::Matrix3cd m_nu = *Dep::m_nu;
      Eigen::Matrix3cd Vnu = *Dep::SeesawI_Vnu;
      Eigen::Matrix3cd Theta = *Dep::SeesawI_Theta;

      vector<double> ml = {sminputs.mE, sminputs.mMu, sminputs.mTau};
      vector<double> mnu = {real(m_nu(0,0)), real(m_nu(1,1)), real(m_nu(2,2)), *Param["M_1"], *Param["M_2"], *Param["M_3"]};
      Eigen::Matrix<complex<double>,3,6> U;

      for(int i=0; i<3; i++)
        for(int j=0; j<3; j++)
        {
          U(i,j) = Vnu(i,j);
          U(i,j+3) = Theta(i,j);
        }

      // Form factors
      int e = 0, mu = 1;
      complex<double> k1r = FormFactors::K1R(mu, e, sminputs, U, mnu);
      complex<double> k2l = FormFactors::K2L(mu, e, sminputs, U, ml, mnu);
      complex<double> k2r = FormFactors::K2R(mu, e, sminputs, U, ml, mnu);

      int u = 0, d =0, s = 1;
      complex<double> CVLLu = FormFactors::CVLL(mu, e, u, u, sminputs, U, ml, mnu);
      complex<double> CVLLd = FormFactors::BVLL(mu, e, d, d, sminputs, U, ml, mnu);
      complex<double> CVLLs = FormFactors::BVLL(mu, e, s, s, sminputs, U, ml, mnu);
      complex<double> CVLRu = FormFactors::CVLR(mu, e, u, u, sminputs, U, ml, mnu);
      complex<double> CVLRd = FormFactors::BVLR(mu, e, d, d, sminputs, U, ml, mnu);
      complex<double> CVLRs = FormFactors::BVLR(mu, e, s, s, sminputs, U, ml, mnu);
      complex<double> CVRLu = FormFactors::CVRL(mu, e, u, u, sminputs, U, ml, mnu);
      complex<double> CVRLd = FormFactors::BVRL(mu, e, d, d, sminputs, U, ml, mnu);
      complex<double> CVRLs = FormFactors::BVRL(mu, e, s, s, sminputs, U, ml, mnu);
      complex<double> CVRRu = FormFactors::CVRR(mu, e, u, u, sminputs, U, ml, mnu);
      complex<double> CVRRd = FormFactors::BVRR(mu, e, d, d, sminputs, U, ml, mnu);
      complex<double> CVRRs = FormFactors::BVRR(mu, e, s, s, sminputs, U, ml, mnu);

      complex<double> CSLLu = FormFactors::CSLL(mu, e, u, u, sminputs, U, ml, mnu, *Param["mH"]);
      complex<double> CSLLd = FormFactors::BSLL(mu, e, d, d, sminputs, U, ml, mnu, *Param["mH"]);
      complex<double> CSLLs = FormFactors::BSLL(mu, e, s, s, sminputs, U, ml, mnu, *Param["mH"]);
      complex<double> CSLRu = FormFactors::CSLL(mu, e, u, u, sminputs, U, ml, mnu, *Param["mH"]);
      complex<double> CSLRd = FormFactors::BSLL(mu, e, d, d, sminputs, U, ml, mnu, *Param["mH"]);
      complex<double> CSLRs = FormFactors::BSLL(mu, e, s, s, sminputs, U, ml, mnu, *Param["mH"]);
      complex<double> CSRLu = FormFactors::CSLL(mu, e, u, u, sminputs, U, ml ,mnu, *Param["mH"]);
      complex<double> CSRLd = FormFactors::BSLL(mu, e, d, d, sminputs, U, ml ,mnu, *Param["mH"]);
      complex<double> CSRLs = FormFactors::BSLL(mu, e, s, s, sminputs, U, ml ,mnu, *Param["mH"]);
      complex<double> CSRRu = FormFactors::CSLL(mu, e, u, u, sminputs, U, ml ,mnu, *Param["mH"]);
      complex<double> CSRRd = FormFactors::BSLL(mu, e, d, d, sminputs, U, ml, mnu, *Param["mH"]);
      complex<double> CSRRs = FormFactors::BSLL(mu, e, s, s, sminputs, U, ml ,mnu, *Param["mH"]);

      double Qu = 2./3.;
      complex<double> gVLu = sqrt(2)/sminputs.GF * (4.*pi / sminputs.alphainv * Qu * (0. - k2r) - 0.5*(CVLLu + CVLRu));
      complex<double> gSLu = -1./(sqrt(2)*sminputs.GF)*(CSLLu + CSLRu);
      complex<double> gVRu = sqrt(2)/sminputs.GF * (4.*pi / sminputs.alphainv * Qu * (k1r - k2l) - 0.5*(CVRRu + CVRLu));
      complex<double> gSRu = -1./(sqrt(2)*sminputs.GF)*(CSRRu + CSRLu);

      double Qd = -1./3.;
      complex<double> gVLd = sqrt(2)/sminputs.GF * (4.*pi / sminputs.alphainv * Qd * (0. - k2r) - 0.5*(CVLLd + CVLRd));
      complex<double> gSLd = -1./(sqrt(2)*sminputs.GF)*(CSLLd + CSLRd);
      complex<double> gVRd = sqrt(2)/sminputs.GF * (4.*pi / sminputs.alphainv * Qd * (k1r - k2l) - 0.5*(CVRRd + CVRLd));
      complex<double> gSRd = -1./(sqrt(2)*sminputs.GF)*(CSRRd + CSRLd);

      double Qs = -1./3.;
      complex<double> gVLs = sqrt(2)/sminputs.GF * (4.*pi / sminputs.alphainv * Qs * (0. - k2r) - 0.5*(CVLLs + CVLRs));
      complex<double> gSLs = -1./(sqrt(2)*sminputs.GF)*(CSLLs + CSLRs);
      complex<double> gVRs = sqrt(2)/sminputs.GF * (4.*pi / sminputs.alphainv * Qs * (k1r - k2l) - 0.5*(CVRRs + CVRLs));
      complex<double> gSRs = -1./(sqrt(2)*sminputs.GF)*(CSRRs + CSRLs);

      double GVup = 2, GVdn = 2, GVdp = 1, GVun = 1, GVsp = 0, GVsn = 0;
      double GSup = 5.1, GSdn = 5.1, GSdp = 4.3, GSun = 4.3, GSsp = 2.5, GSsn = 2.5;
      complex<double> g0SL = 0.5*(gSLu*(GSup + GSun) + gSLd*(GSdp + GSdn) + gSLs*(GSsp + GSsn));
      complex<double> g0SR = 0.5*(gSRu*(GSup + GSun) + gSRd*(GSdp + GSdn) + gSRs*(GSsp + GSsn));
      complex<double> g0VL = 0.5*(gVLu*(GVup + GVun) + gVLd*(GVdp + GVdn) + gVLs*(GVsp + GVsn));
      complex<double> g0VR = 0.5*(gVRu*(GVup + GVun) + gVRd*(GVdp + GVdn) + gVRs*(GVsp + GVsn));
      complex<double> g1SL = 0.5*(gSLu*(GSup - GSun) + gSLd*(GSdp - GSdn) + gSLs*(GSsp - GSsn));
      complex<double> g1SR = 0.5*(gSRu*(GSup - GSun) + gSRd*(GSdp - GSdn) + gSRs*(GSsp - GSsn));
      complex<double> g1VL = 0.5*(gVLu*(GVup - GVun) + gVLd*(GVdp - GVdn) + gVLs*(GVsp - GVsn));
      complex<double> g1VR = 0.5*(gVRu*(GVup - GVun) + gVRd*(GVdp - GVdn) + gVRs*(GVsp - GVsn));


      // Parameters for Pb, from Table 1 in 1209.2679 for Pb
      double Z = 82, N = 126;
      double Zeff = 34., Fp = 0.15;
      double hbar = 6.582119514e-25; // GeV * s
      double GammaCapt = 13.45e6 * hbar;

      result = (pow(sminputs.GF,2)*pow(sminputs.mMu,5)*pow(Zeff,4)*pow(Fp,2)) / (8.*pow(pi,4)*pow(sminputs.alphainv,3)*Z*GammaCapt) * (norm((Z+N)*(g0VL + g0SL) + (Z-N)*(g1VL + g1SL)) + norm((Z+N)*(g0VR + g0SR) + (Z-N)*(g1VR + g1SR)));

    }

    /// Likelihood for l -> l gamma processes
    void l2lgamma_likelihood(double &result)
    {
      using namespace Pipes::l2lgamma_likelihood;

      static bool first = true;
      static boost::numeric::ublas::matrix<double> cov_exp, value_exp;
      static double th_err[3];
      double theory[3];

      // Read and calculate things based on the observed data only the first time through, as none of it depends on the model parameters.
      if (first)
      {
        // Read in experimental measuremens
        Flav_reader fread(GAMBIT_DIR  "/FlavBit/data");
        fread.debug_mode(flav_debug);

        // mu -> e gamma
        fread.read_yaml_measurement("flav_data.yaml", "BR_muegamma");
        // tau -> e gamma
        fread.read_yaml_measurement("flav_data.yaml", "BR_tauegamma");
        // tau -> mu gamma
        fread.read_yaml_measurement("flav_data.yaml", "BR_taumugamma");

        fread.initialise_matrices();
        cov_exp=fread.get_exp_cov();
        value_exp=fread.get_exp_value();

        for (int i = 0; i < 3; ++i)
          th_err[i] = fread.get_th_err()(i,0).first;

        // Init over.
        first = false;
      }

     theory[0] = *Dep::muegamma;
     if(flav_debug) cout << "mu- -> e- gamma = " << theory[0] << endl;
     theory[1] = *Dep::tauegamma;
     if(flav_debug) cout << "tau- -> e- gamma = " << theory[1] << endl;
     theory[2] = *Dep::taumugamma;
     if(flav_debug) cout << "tau- -> mu- gamma = " << theory[2] << endl;

     result = 0;
     for (int i = 0; i < 3; ++i)
       result += Stats::gaussian_upper_limit(theory[i], value_exp(i,0), th_err[i], sqrt(cov_exp(i,i)), false);

    }

    /// Likelihood for l -> l l l processes
    void l2lll_likelihood(double &result)
    {
      using namespace Pipes::l2lll_likelihood;

      static bool first = true;
      static boost::numeric::ublas::matrix<double> cov_exp, value_exp;
      static double th_err[7];
      double theory[7];


      // Read and calculate things based on the observed data only the first time through, as none of it depends on the model parameters.
      if (first)
      {
        // Read in experimental measuremens
        Flav_reader fread(GAMBIT_DIR  "/FlavBit/data");
        fread.debug_mode(flav_debug);

        // mu- -> e- e- e+
        fread.read_yaml_measurement("flav_data.yaml", "BR_mueee");
        // tau- -> e- e- e+
        fread.read_yaml_measurement("flav_data.yaml", "BR_taueee");
        // tau- -> mu- mu- mu+
        fread.read_yaml_measurement("flav_data.yaml", "BR_taumumumu");
        // tau- -> mu- e- e+
        fread.read_yaml_measurement("flav_data.yaml", "BR_taumuee");
        // tau- -> e- e- mu+
        fread.read_yaml_measurement("flav_data.yaml", "BR_taueemu");
        // tau- -> e- mu- mu+
        fread.read_yaml_measurement("flav_data.yaml", "BR_tauemumu");
        // tau- -> mu- mu- e+
        fread.read_yaml_measurement("flav_data.yaml", "BR_taumumue");

        fread.initialise_matrices();
        cov_exp=fread.get_exp_cov();
        value_exp=fread.get_exp_value();

        for (int i = 0; i < 7; ++i)
          th_err[i] = fread.get_th_err()(i,0).first;

        // Init over.
        first = false;
      }

     theory[0] = *Dep::mueee;
     if(flav_debug) cout << "mu-  -> e-  e-  e+  = " << theory[0] << endl;
     theory[1] = *Dep::taueee;
     if(flav_debug) cout << "tau- -> e-  e-  e+  = " << theory[1] << endl;
     theory[2] = *Dep::taumumumu;
     if(flav_debug) cout << "tau- -> mu- mu- mu+ = " << theory[2] << endl;
     theory[3] = *Dep::taumuee;
     if(flav_debug) cout << "tau- -> mu- e-  e-  = " << theory[3] << endl;
     theory[4] = *Dep::taueemu;
     if(flav_debug) cout << "tau- -> e-  e-  mu+ = " << theory[4] << endl;
     theory[5] = *Dep::tauemumu;
     if(flav_debug) cout << "tau- -> e-  mu- mu+ = " << theory[5] << endl;
     theory[6] = *Dep::taumumue;
     if(flav_debug) cout << "tau- -> mu- mu- e+  = " << theory[6] << endl;

     result = 0;
     for (int i = 0; i < 7; ++i)
       result += Stats::gaussian_upper_limit(theory[i], value_exp(i,0), th_err[i], sqrt(cov_exp(i,i)), false);

    }

    /// Likelihood for mu - e conversion in nucleii
    void mu2e_likelihood(double &result)
    {
      using namespace Pipes::mu2e_likelihood;

      static bool first = true;
      static boost::numeric::ublas::matrix<double> cov_exp, value_exp;
      static double th_err[2];
      double theory[2];


      // Read and calculate things based on the observed data only the first time through, as none of it depends on the model parameters.
      if (first)
      {
        // Read in experimental measuremens
        Flav_reader fread(GAMBIT_DIR  "/FlavBit/data");
        fread.debug_mode(flav_debug);

        // mu - e (Ti)
        fread.read_yaml_measurement("flav_data.yaml", "R_mueTi");
        // mu - e (Pb)
        fread.read_yaml_measurement("flav_data.yaml", "R_muePb");

        fread.initialise_matrices();
        cov_exp=fread.get_exp_cov();
        value_exp=fread.get_exp_value();

        for (int i = 0; i < 2; ++i)
          th_err[i] = fread.get_th_err()(i,0).first;

        // Init over.
        first = false;
      }

      theory[0] = *Dep::mueTi;
      if(flav_debug) cout << "mu - e (Ti) = " << theory[0] << endl;
      theory[1] = *Dep::muePb;
      if(flav_debug) cout << "mu - e (Pb) = " << theory[1] << endl;

      result = 0;
      for (int i = 0; i < 2; ++i)
        result += Stats::gaussian_upper_limit(theory[i], value_exp(i,0), th_err[i], sqrt(cov_exp(i,i)), false);

    }

    /// Measurements for LUV in b->sll
    void LUV_measurements(predictions_measurements_covariances &pmc)
    {
      using namespace Pipes::LUV_measurements;
      static bool first = true;

      static double theory_RKstar_0045_11_err, theory_RKstar_11_60_err, theory_RK_err;
      if (flav_debug) cout<<"Starting LUV_measurements"<<endl;

      // Read and calculate things based on the observed data only the first time through, as none of it depends on the model parameters.
      if (first)
        {
          pmc.LL_name="LUV_likelihood";

          Flav_reader fread(GAMBIT_DIR  "/FlavBit/data");
          fread.debug_mode(flav_debug);

          if (flav_debug) cout<<"Initiated Flav reader in LUV_measurements"<<endl;
          fread.read_yaml_measurement("flav_data.yaml", "RKstar_0045_11");
          fread.read_yaml_measurement("flav_data.yaml", "RKstar_11_60");
          fread.read_yaml_measurement("flav_data.yaml", "RK");

          if (flav_debug) cout<<"Finished reading LUV data"<<endl;

          fread.initialise_matrices();

          theory_RKstar_0045_11_err = fread.get_th_err()(0,0).first;
          theory_RKstar_11_60_err = fread.get_th_err()(1,0).first;
          theory_RK_err = fread.get_th_err()(2,0).first;

          pmc.value_exp=fread.get_exp_value();
          pmc.cov_exp=fread.get_exp_cov();

          pmc.value_th.resize(3,1);
          pmc.cov_th.resize(3,3);

          pmc.dim=3;

          // Init over and out.
          first = false;
        }

      // Get theory prediction
      pmc.value_th(0,0)=*Dep::RKstar_0045_11;
      pmc.value_th(1,0)=*Dep::RKstar_11_60;
      pmc.value_th(2,0)=*Dep::RK;

      // Compute error on theory prediction and populate the covariance matrix
      pmc.cov_th(0,0)=theory_RKstar_0045_11_err;
      pmc.cov_th(0,1)=0.;
      pmc.cov_th(0,2)=0.;
      pmc.cov_th(1,0)=0.;
      pmc.cov_th(1,1)=theory_RKstar_11_60_err;
      pmc.cov_th(1,2)=0.;
      pmc.cov_th(2,0)=0.;
      pmc.cov_th(2,1)=0.;
      pmc.cov_th(2,2)=theory_RK_err;



      // Save the differences between theory and experiment
      pmc.diff.clear();
      for (int i=0;i<3;++i)
        {
          pmc.diff.push_back(pmc.value_exp(i,0)-pmc.value_th(i,0));
        }

      if (flav_debug) cout<<"Finished LUV_measurements"<<endl;


    }
    /// Likelihood  for LUV in b->sll
    void LUV_likelihood(double &result)
    {
      using namespace Pipes::LUV_likelihood;

      if (flav_debug) cout<<"Starting LUV_likelihood"<<endl;

      predictions_measurements_covariances pmc = *Dep::LUV_M;

      boost::numeric::ublas::matrix<double> cov=pmc.cov_exp;

      // adding theory and experimental covariance
      cov+=pmc.cov_th;

      //calculating a diff
      vector<double> diff;
      diff=pmc.diff;

      boost::numeric::ublas::matrix<double> cov_inv(pmc.dim, pmc.dim);
      InvertMatrix(cov, cov_inv);

      double Chi2=0;
      for (int i=0; i < pmc.dim; ++i)
        {
          for (int j=0; j<pmc.dim; ++j)
            {
              Chi2+= diff[i] * cov_inv(i,j)*diff[j];
            }
        }

      result=-0.5*Chi2;

      if (flav_debug) cout<<"Finished LUV_likelihood"<<endl;

      if (flav_debug_LL) cout<<"Likelihood result LUV_likelihood  : "<< result<<endl;

    }

    /// Br Bs->mumu decays for the untagged case (CP-averaged)
    void Flavio_test(double &result)
    {
      using namespace Pipes::Flavio_test;
      if (flav_debug) cout<<"Starting Flavio_test"<<endl;

      result=BEreq::sm_prediction_CONV("BR(Bs->mumu)");
      std::cout<<"Flavio result: "<<result<<std::endl;
    }



    /// HEPLike LogLikelihood B -> tau nu
    void HEPLike_B2TauNuLogLikelihood(double &result)
    {
      using namespace Pipes::HEPLike_B2TauNuLogLikelihood;
      static const std::string inputfile = path_to_latest_heplike_data() + "/data/PDG/Semileptonic/B2TauNu.yaml";
      static HepLike_default::HL_Gaussian gaussian(inputfile);
      static bool first = true;
      if (first)
      {
        std::cout << "Debug: Reading HepLike data file: " << inputfile << endl;
        gaussian.Read();
        first = false;
      }
      const double theory = *Dep::Btaunu;
      result = gaussian.GetLogLikelihood(theory /* , theory_error */);
      if (flav_debug) std::cout << "hepLikeB2TauNuLogLikelihood result: " << result << std::endl;
    }

    /// HEPLike LogLikelihood RD RDstar
    void HEPLike_RDRDstarLogLikelihood(double& result)
    {
      using namespace Pipes::HEPLike_RDRDstarLogLikelihood;
      static const std::string inputfile = path_to_latest_heplike_data() + "/data/HFLAV_18/Semileptonic/RD_RDstar.yaml";
      static HepLike_default::HL_nDimGaussian nDimGaussian(inputfile);
      static bool first = true;
      if (first)
      {
        std::cout << "Debug: Reading HepLike data file: " << inputfile << endl;
        nDimGaussian.Read();
        first = false;
      }
      const std::vector<double> theory{*Dep::RD, *Dep::RDstar};
      result = nDimGaussian.GetLogLikelihood(theory /* , theory_covariance */);
      if (flav_debug) std::cout << "hepLikeRDRDstarLogLikelihood result: " << result << std::endl;
    }

    /// HEPLike LogLikelihood b -> s gamma
    void HEPLike_B2SGammaLogLikelihood(double &result)
    {
      using namespace Pipes::HEPLike_B2SGammaLogLikelihood;
      static const std::string inputfile = path_to_latest_heplike_data() + "/data/HFLAV_18/RD/b2sgamma.yaml";
      static HepLike_default::HL_Gaussian gaussian(inputfile);
      static bool first = true;
      if (first)
      {
        std::cout << "Debug: Reading HepLike data file: " << inputfile << endl;
        gaussian.Read();
        first = false;
      }

      static const std::string observable{"BR_BXsgamma"};

      auto SI_theory = *Dep::SuperIso_obs_values;
      auto SI_theory_covariance = *Dep::SuperIso_theory_covariance;

      result = gaussian.GetLogLikelihood(
              SI_theory[observable],
              SI_theory_covariance[observable][observable]
              );
      if (flav_debug) std::cout << "hepLikeB2SGammaLogLikelihood result: " << result << std::endl;
    }

    /// HEPLike LogLikelihood B -> ll
    void HEPLike_B2mumuLogLikelihood(double &result)
    {
      using namespace Pipes::HEPLike_B2mumuLogLikelihood;
      static const std::string inputfile_LHCb = path_to_latest_heplike_data() + "/data/LHCb/RD/B2MuMu/CERN-EP-2017-100.yaml";
      static const std::string inputfile_CMS = path_to_latest_heplike_data() + "/data/CMS/RD/B2MuMu/CERN-EP-2017-100.yaml"; // THIS NEEDS TO BE IMPLEMENTED

      static HepLike_default::HL_nDimLikelihood nDimLikelihood(inputfile_LHCb);
      //static HepLike_default::HL_nDimLikelihood nDimLikelihood(inputfile_CMS);

      static bool first = true;
      if (first)
      {
        std::cout << "Debug: Reading HepLike data file: " << inputfile_LHCb << endl;
        nDimLikelihood.Read();

        first = false;
      }

      static const std::array<std::string, 2> observables{
        "BRuntag_Bsmumu",
        "BR_Bdmumu"
      };

      auto SI_theory = *Dep::SuperIso_obs_values;
      auto SI_theory_covariance = *Dep::SuperIso_theory_covariance;

      // C++14 allows auto instead of decltype(observables0p1_0p98)
      auto get_obs_theory = [SI_theory](decltype(observables)& observables){
          std::vector<double> obs_theory;
          obs_theory.reserve(observables.size());
          for (unsigned int i = 0; i < observables.size(); ++i) {
            obs_theory.push_back(SI_theory.at(observables[i]));
          }
          return obs_theory;
      };

      auto get_obs_covariance = [SI_theory_covariance](decltype(observables)& observables){
          boost::numeric::ublas::matrix<double> obs_covariance(observables.size(), observables.size());
          for (unsigned int i = 0; i < observables.size(); ++i) {
            for (unsigned int j = 0; j < observables.size(); ++j) {
              obs_covariance(i, j) = SI_theory_covariance.at(observables[i]).at(observables[j]);
            }
          }
          return obs_covariance;
      };

      result = nDimLikelihood.GetLogLikelihood(
              get_obs_theory(observables)
              /* nDimLikelihood does not support theory errors */
              );

      if (flav_debug) std::cout << "%s result: " << result << std::endl;
    }

    /// HEPLike LogLikelihood B -> K* mu mu Angular
    void HEPLike_B2KstarmumuAng_LogLikelihood(double &result)
    {
      using namespace Pipes::HEPLike_B2KstarmumuAng_LogLikelihood;

      static const std::string inputfile_q2_0p1_1p1 = path_to_latest_heplike_data() + "/data/LHCb/RD/Bd2KstarMuMu_Angular/PH-EP-2015-314_q2_0.1_0.98.yaml";
      static const std::string inputfile_q2_1p1_2p5 = path_to_latest_heplike_data() + "/data/LHCb/RD/Bd2KstarMuMu_Angular/PH-EP-2015-314_q2_1.1_2.5.yaml";
      static const std::string inputfile_q2_2p5_4 = path_to_latest_heplike_data() + "/data/LHCb/RD/Bd2KstarMuMu_Angular/PH-EP-2015-314_q2_2.5_4.0.yaml";
      static const std::string inputfile_q2_4_6 = path_to_latest_heplike_data() + "/data/LHCb/RD/Bd2KstarMuMu_Angular/PH-EP-2015-314_q2_4.0_6.0.yaml";
      static const std::string inputfile_q2_6_8 = path_to_latest_heplike_data() + "/data/LHCb/RD/Bd2KstarMuMu_Angular/PH-EP-2015-314_q2_6.0_8.0.yaml";
      static const std::string inputfile_q2_15_19 = path_to_latest_heplike_data() + "/data/LHCb/RD/Bd2KstarMuMu_Angular/PH-EP-2015-314_q2_15.0_19.yaml";
      static HepLike_default::HL_nDimBifurGaussian nDimBifurGaussian_0(inputfile_q2_0p1_1p1);
      static HepLike_default::HL_nDimBifurGaussian nDimBifurGaussian_1(inputfile_q2_1p1_2p5);
      static HepLike_default::HL_nDimBifurGaussian nDimBifurGaussian_2(inputfile_q2_2p5_4);
      static HepLike_default::HL_nDimBifurGaussian nDimBifurGaussian_3(inputfile_q2_4_6);
      static HepLike_default::HL_nDimBifurGaussian nDimBifurGaussian_4(inputfile_q2_6_8);
      static HepLike_default::HL_nDimBifurGaussian nDimBifurGaussian_5(inputfile_q2_15_19);

      static bool first = true;
      if (first)
      {
        std::cout << "Debug: Reading HepLike data file: " << inputfile_q2_0p1_1p1 << endl;
        nDimBifurGaussian_0.Read();
        std::cout << "Debug: Reading HepLike data file: " << inputfile_q2_1p1_2p5 << endl;
        nDimBifurGaussian_1.Read();
        std::cout << "Debug: Reading HepLike data file: " << inputfile_q2_2p5_4 << endl;
        nDimBifurGaussian_2.Read();
        std::cout << "Debug: Reading HepLike data file: " << inputfile_q2_4_6 << endl;
        nDimBifurGaussian_3.Read();
        std::cout << "Debug: Reading HepLike data file: " << inputfile_q2_6_8 << endl;
        nDimBifurGaussian_4.Read();
        std::cout << "Debug: Reading HepLike data file: " << inputfile_q2_15_19 << endl;
        nDimBifurGaussian_5.Read();

        first = false;
      }

      // Ordering of observables defined by HEPLike
      static const std::array<std::string, 8> observables0p1_0p98{
        "FL_B0Kstar0mumu_0.1_0.98",
        "S3_B0Kstar0mumu_0.1_0.98",
        "S4_B0Kstar0mumu_0.1_0.98",
        "S5_B0Kstar0mumu_0.1_0.98",
        "AFB_B0Kstar0mumu_0.1_0.98",
        "S7_B0Kstar0mumu_0.1_0.98",
        "S8_B0Kstar0mumu_0.1_0.98",
        "S9_B0Kstar0mumu_0.1_0.98",
      };
      static const std::array<std::string, 8> observables1p1_2p5{
        "FL_B0Kstar0mumu_1.1_2.5",
        "S3_B0Kstar0mumu_1.1_2.5",
        "S4_B0Kstar0mumu_1.1_2.5",
        "S5_B0Kstar0mumu_1.1_2.5",
        "AFB_B0Kstar0mumu_1.1_2.5",
        "S7_B0Kstar0mumu_1.1_2.5",
        "S8_B0Kstar0mumu_1.1_2.5",
        "S9_B0Kstar0mumu_1.1_2.5",
      };
      static const std::array<std::string, 8> observables2p5_4{
        "FL_B0Kstar0mumu_2.5_4",
        "S3_B0Kstar0mumu_2.5_4",
        "S4_B0Kstar0mumu_2.5_4",
        "S5_B0Kstar0mumu_2.5_4",
        "AFB_B0Kstar0mumu_2.5_4",
        "S7_B0Kstar0mumu_2.5_4",
        "S8_B0Kstar0mumu_2.5_4",
        "S9_B0Kstar0mumu_2.5_4",
      };
      static const std::array<std::string, 8> observables4_6{
        "FL_B0Kstar0mumu_4_6",
        "S3_B0Kstar0mumu_4_6",
        "S4_B0Kstar0mumu_4_6",
        "S5_B0Kstar0mumu_4_6",
        "AFB_B0Kstar0mumu_4_6",
        "S7_B0Kstar0mumu_4_6",
        "S8_B0Kstar0mumu_4_6",
        "S9_B0Kstar0mumu_4_6",
      };
      static const std::array<std::string, 8> observables6_8{
        "FL_B0Kstar0mumu_6_8",
        "S3_B0Kstar0mumu_6_8",
        "S4_B0Kstar0mumu_6_8",
        "S5_B0Kstar0mumu_6_8",
        "AFB_B0Kstar0mumu_6_8",
        "S7_B0Kstar0mumu_6_8",
        "S8_B0Kstar0mumu_6_8",
        "S9_B0Kstar0mumu_6_8",
      };
      static const std::array<std::string, 8> observables15_19{
        "FL_B0Kstar0mumu_15_19",
        "S3_B0Kstar0mumu_15_19",
        "S4_B0Kstar0mumu_15_19",
        "S5_B0Kstar0mumu_15_19",
        "AFB_B0Kstar0mumu_15_19",
        "S7_B0Kstar0mumu_15_19",
        "S8_B0Kstar0mumu_15_19",
        "S9_B0Kstar0mumu_15_19",
      };

      auto SI_theory = *Dep::SuperIso_obs_values;
      auto SI_theory_covariance = *Dep::SuperIso_theory_covariance;
      
      // C++14 allows auto instead of decltype(observables0p1_0p98)
      auto get_obs_theory = [SI_theory](decltype(observables0p1_0p98)& observables){
        std::vector<double> obs_theory;
        obs_theory.reserve(observables.size());
        for (unsigned int i = 0; i < observables.size(); ++i) {
          obs_theory.push_back(SI_theory.at(observables[i]));
        }
        return obs_theory;
      };

      auto get_obs_covariance = [SI_theory_covariance](decltype(observables0p1_0p98)& observables){
        boost::numeric::ublas::matrix<double> obs_covariance(observables.size(), observables.size());
        for (unsigned int i = 0; i < observables.size(); ++i) {
          for (unsigned int j = 0; j < observables.size(); ++j) {
            obs_covariance(i, j) = SI_theory_covariance.at(observables[i]).at(observables[j]);
          }
        }
        return obs_covariance;
      };
      result = 0;
      result += nDimBifurGaussian_0.GetLogLikelihood(get_obs_theory(observables0p1_0p98), get_obs_covariance(observables0p1_0p98));
      result += nDimBifurGaussian_1.GetLogLikelihood(get_obs_theory(observables1p1_2p5), get_obs_covariance(observables1p1_2p5));
      result += nDimBifurGaussian_2.GetLogLikelihood(get_obs_theory(observables2p5_4), get_obs_covariance(observables2p5_4));
      result += nDimBifurGaussian_3.GetLogLikelihood(get_obs_theory(observables4_6), get_obs_covariance(observables4_6));
      result += nDimBifurGaussian_4.GetLogLikelihood(get_obs_theory(observables6_8), get_obs_covariance(observables6_8));
      result += nDimBifurGaussian_5.GetLogLikelihood(get_obs_theory(observables15_19), get_obs_covariance(observables15_19));

      if (flav_debug) std::cout << "%s result: " << result << std::endl;
    }

    /// HEPLike LogLikelihood B -> K* mu mu Br
    void HEPLike_B2KstarmumuBr_LogLikelihood(double &result)
    {

      using namespace Pipes::HEPLike_B2KstarmumuBr_LogLikelihood;

      static const std::string inputfile_q2_0p1_1p1 = path_to_latest_heplike_data() + "/data/LHCb/RD/Bd2KstarMuMu_Br/CERN-EP-2016-141_q2_0.1_0.98.yaml";
      static const std::string inputfile_q2_1p1_2p5 = path_to_latest_heplike_data() + "/data/LHCb/RD/Bd2KstarMuMu_Br/CERN-EP-2016-141_q2_1.1_2.5.yaml";
      static const std::string inputfile_q2_2p5_4 = path_to_latest_heplike_data() + "/data/LHCb/RD/Bd2KstarMuMu_Br/CERN-EP-2016-141_q2_2.5_4.yaml";
      static const std::string inputfile_q2_4_6 = path_to_latest_heplike_data() + "/data/LHCb/RD/Bd2KstarMuMu_Br/CERN-EP-2016-141_q2_4_6.yaml";
      static const std::string inputfile_q2_6_8 = path_to_latest_heplike_data() + "/data/LHCb/RD/Bd2KstarMuMu_Br/CERN-EP-2016-141_q2_6_8.yaml";
      static const std::string inputfile_q2_15_19 = path_to_latest_heplike_data() + "/data/LHCb/RD/Bd2KstarMuMu_Br/CERN-EP-2016-141_q2_15_19.yaml";
      static HepLike_default::HL_BifurGaussian BifurGaussian_0(inputfile_q2_0p1_1p1);
      static HepLike_default::HL_BifurGaussian BifurGaussian_1(inputfile_q2_1p1_2p5);
      static HepLike_default::HL_BifurGaussian BifurGaussian_2(inputfile_q2_2p5_4);
      static HepLike_default::HL_BifurGaussian BifurGaussian_3(inputfile_q2_4_6);
      static HepLike_default::HL_BifurGaussian BifurGaussian_4(inputfile_q2_6_8);
      static HepLike_default::HL_BifurGaussian BifurGaussian_5(inputfile_q2_15_19);

      static bool first = true;
      if (first)
      {
        std::cout << "Debug: Reading HepLike data file: " << inputfile_q2_0p1_1p1 << endl;
        BifurGaussian_0.Read();
        std::cout << "Debug: Reading HepLike data file: " << inputfile_q2_1p1_2p5 << endl;
        BifurGaussian_1.Read();
        std::cout << "Debug: Reading HepLike data file: " << inputfile_q2_2p5_4 << endl;
        BifurGaussian_2.Read();
        std::cout << "Debug: Reading HepLike data file: " << inputfile_q2_4_6 << endl;
        BifurGaussian_3.Read();
        std::cout << "Debug: Reading HepLike data file: " << inputfile_q2_6_8 << endl;
        BifurGaussian_4.Read();
        std::cout << "Debug: Reading HepLike data file: " << inputfile_q2_15_19 << endl;
        BifurGaussian_5.Read();

        first = false;
      }

      // Ordering of observables defined by HEPLike
      // Nota bene: Although the variables are called dGamma, these functions actually return the differential BR.
      //            This holds true for SuperIso 4.1, could change in future versions though.
      static const std::array<std::string, 6> observables{
        "dGamma/dq2_B0Kstar0mumu_0.1_0.98",
        "dGamma/dq2_B0Kstar0mumu_1.1_2.5",
        "dGamma/dq2_B0Kstar0mumu_2.5_4",
        "dGamma/dq2_B0Kstar0mumu_4_6",
        "dGamma/dq2_B0Kstar0mumu_6_8",
        "dGamma/dq2_B0Kstar0mumu_15_19",
      };

      auto SI_theory = *Dep::SuperIso_obs_values;
      auto SI_theory_covariance = *Dep::SuperIso_theory_covariance;

      result = 0;
      result += BifurGaussian_0.GetLogLikelihood(SI_theory[observables[0]], SI_theory_covariance[observables[0]][observables[0]]);
      result += BifurGaussian_1.GetLogLikelihood(SI_theory[observables[1]], SI_theory_covariance[observables[1]][observables[1]]);
      result += BifurGaussian_2.GetLogLikelihood(SI_theory[observables[2]], SI_theory_covariance[observables[2]][observables[2]]);
      result += BifurGaussian_3.GetLogLikelihood(SI_theory[observables[3]], SI_theory_covariance[observables[3]][observables[3]]);
      result += BifurGaussian_4.GetLogLikelihood(SI_theory[observables[4]], SI_theory_covariance[observables[4]][observables[4]]);
      result += BifurGaussian_5.GetLogLikelihood(SI_theory[observables[5]], SI_theory_covariance[observables[5]][observables[5]]);

      if (flav_debug) std::cout << "%s result: " << result << std::endl;
    }

    void HEPLike_Bs2PhimumuBr_LogLikelihood(double &result)
    {
      using namespace Pipes::HEPLike_Bs2PhimumuBr_LogLikelihood;

      static const std::string inputfile_0 = path_to_latest_heplike_data() + "/data/LHCb/RD/Bs2PhiMuMu_Br/CERN-PH-EP-2015-145_1_6.yaml";
      static const std::string inputfile_1 = path_to_latest_heplike_data() + "/data/LHCb/RD/Bs2PhiMuMu_Br/CERN-PH-EP-2015-145_15_19.yaml";
      static HepLike_default::HL_BifurGaussian bifurGaussian_0(inputfile_0);
      static HepLike_default::HL_BifurGaussian bifurGaussian_1(inputfile_1);

      static bool first = true;
      if (first)
      {
        std::cout << "Debug: Reading HepLike data file: " << inputfile_0 << endl;
        bifurGaussian_0.Read();
        std::cout << "Debug: Reading HepLike data file: " << inputfile_1 << endl;
        bifurGaussian_1.Read();

        first = false;
      }

      // Ordering of observables defined by HEPLike
      // Nota bene: Although the variables are called dGamma, these functions actually return the differential BR.
      //            This holds true for SuperIso 4.1, could change in future versions though.
      static const std::array<std::string, 2> observables{
              "dGamma/dq2_Bsphimumu_1_6",
              "dGamma/dq2_Bsphimumu_15_19",
      };

      auto SI_theory = *Dep::SuperIso_obs_values;
      auto SI_theory_covariance = *Dep::SuperIso_theory_covariance;

      result = 0;
      result += bifurGaussian_0.GetLogLikelihood(SI_theory[observables[0]], SI_theory_covariance[observables[0]][observables[0]]);
      result += bifurGaussian_1.GetLogLikelihood(SI_theory[observables[1]], SI_theory_covariance[observables[1]][observables[1]]);

      if (flav_debug) std::cout << "%s result: " << result << std::endl;
    }
  }
}
