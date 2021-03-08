//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  Implementation of MDMSM
///  DarkBit routines.
///
///  Authors (add name and date if you modify):    
///       *** Automatically created by GUM ***     
///                                                
///  \author The GAMBIT Collaboration             
///  \date 04:53PM on March 08, 2021
///                                                
///  ********************************************* 

#include "gambit/Elements/gambit_module_headers.hpp"
#include "gambit/DarkBit/DarkBit_rollcall.hpp"
#include "gambit/Utils/ascii_table_reader.hpp"
#include "boost/make_shared.hpp"
#include "gambit/DarkBit/DarkBit_utils.hpp"

namespace Gambit
{
  namespace DarkBit
  {
    class MDMSM
    {
      public:
      /// Initialize MDMSM object (branching ratios etc)
      MDMSM() {};
      ~MDMSM() {};
      
      // Annihilation cross-section. sigmav is a pointer to a CalcHEP backend function.
      double sv(str channel, DecayTable& tbl, double (*sigmav)(str&, std::vector<str>&, std::vector<str>&, double&, const DecayTable&), double v_rel)
      {
        /// Returns sigma*v for a given channel.
        double GeV2tocm3s1 = gev2cm2*s2cm;
        
        // CalcHEP args
        str model = "MDMSM"; // CalcHEP model name
        std::vector<str> in = {"~chi", "~chi"}; // In states: DM+DMbar
        std::vector<str> out; // Out states
        if (channel == "dbar_1, d_1") out = {"d~", "d"};
        if (channel == "ubar_2, u_2") out = {"c~", "c"};
        if (channel == "g, g") out = {"g", "g"};
        if (channel == "dbar_3, d_3") out = {"b~", "b"};
        if (channel == "e+_1, e-_1") out = {"e+", "e-"};
        if (channel == "ubar_1, u_1") out = {"u~", "u"};
        if (channel == "e+_2, e-_2") out = {"mu+", "mu-"};
        if (channel == "Y, Y") out = {"Y", "Y"};
        if (channel == "ubar_3, u_3") out = {"t~", "t"};
        if (channel == "dbar_2, d_2") out = {"s~", "s"};
        
        // Check the channel has been filled
        if (out.size() > 1) return sigmav(model, in, out, v_rel, tbl)*GeV2tocm3s1;
        else return 0;
      }
      
      
    };
    
    void TH_ProcessCatalog_MDMSM(TH_ProcessCatalog &result)
    {
      using namespace Pipes::TH_ProcessCatalog_MDMSM;
      using std::vector;
      using std::string;
      
      // Initialize empty catalog, main annihilation process
      TH_ProcessCatalog catalog;
      TH_Process process_ann("~chi", "~chi");
      
      // Explicitly state that DM is self-conjugate
      process_ann.isSelfConj = true;
      
      
      // Import particle masses 
      
      // Convenience macros
      #define getSMmass(Name, spinX2) catalog.particleProperties.insert(std::pair<string, TH_ParticleProperty> (Name, TH_ParticleProperty(SM.get(Par::Pole_Mass,Name), spinX2)));
      #define addParticle(Name, Mass, spinX2) catalog.particleProperties.insert(std::pair<string, TH_ParticleProperty> (Name, TH_ParticleProperty(Mass, spinX2)));
      
      // Import Spectrum objects
      const Spectrum& spec = *Dep::MDMSM_spectrum;
      const SubSpectrum& SM = spec.get_LE();
      const SMInputs& SMI   = spec.get_SMInputs();
      
      // Get SM pole masses
      getSMmass("e-_1",     1)
      getSMmass("e+_1",     1)
      getSMmass("e-_2",     1)
      getSMmass("e+_2",     1)
      getSMmass("e-_3",     1)
      getSMmass("e+_3",     1)
      getSMmass("Z0",       2)
      getSMmass("W+",       2)
      getSMmass("W-",       2)
      getSMmass("g",        2)
      getSMmass("gamma",    2)
      getSMmass("u_3",      1)
      getSMmass("ubar_3",   1)
      getSMmass("d_3",      1)
      getSMmass("dbar_3",   1)
      
      // Pole masses not available for the light quarks.
      addParticle("u_1"   , SMI.mU,  1) // mu(2 GeV)^MS-bar
      addParticle("ubar_1", SMI.mU,  1) // mu(2 GeV)^MS-bar
      addParticle("d_1"   , SMI.mD,  1) // md(2 GeV)^MS-bar
      addParticle("dbar_1", SMI.mD,  1) // md(2 GeV)^MS-bar
      addParticle("u_2"   , SMI.mCmC,1) // mc(mc)^MS-bar
      addParticle("ubar_2", SMI.mCmC,1) // mc(mc)^MS-bar
      addParticle("d_2"   , SMI.mS,  1) // ms(2 GeV)^MS-bar
      addParticle("dbar_2", SMI.mS,  1) // ms(2 GeV)^MS-bar
      
      // Masses for neutrino flavour eigenstates. Set to zero.
      // (presently not required)
      addParticle("nu_e",     0.0, 1)
      addParticle("nubar_e",  0.0, 1)
      addParticle("nu_mu",    0.0, 1)
      addParticle("nubar_mu", 0.0, 1)
      addParticle("nu_tau",   0.0, 1)
      addParticle("nubar_tau",0.0, 1)
      
      // Meson masses
      addParticle("pi0",   meson_masses.pi0,       0)
      addParticle("pi+",   meson_masses.pi_plus,   0)
      addParticle("pi-",   meson_masses.pi_minus,  0)
      addParticle("eta",   meson_masses.eta,       0)
      addParticle("rho0",  meson_masses.rho0,      1)
      addParticle("rho+",  meson_masses.rho_plus,  1)
      addParticle("rho-",  meson_masses.rho_minus, 1)
      addParticle("omega", meson_masses.omega,     1)
      
      // MDMSM-specific masses
      double mchi = spec.get(Par::Pole_Mass, "~chi");
      addParticle("~chi", mchi, 1);
      addParticle("h0_1", spec.get(Par::Pole_Mass, "h0_1"), 0);
      addParticle("Y", spec.get(Par::Pole_Mass, "Y"), 0);
      
      // Get rid of convenience macros
      #undef getSMmass
      #undef addParticle
      
      // Import decay table from DecayBit
      DecayTable tbl = *Dep::decay_rates;
      
      // Set of imported decays
      std::set<string> importedDecays;
      
      // Minimum branching ratio to include
      double minBranching = runOptions->getValueOrDef<double>(0.0, "ProcessCatalog_MinBranching");
      
      // Import relevant decays
      using DarkBit_utils::ImportDecays;
      
      auto excludeDecays = daFunk::vec<std::string>("Z0", "W+", "W-", "e+_3", "e-_3", "e+_2", "e-_2");
      
      ImportDecays("h0_1", catalog, importedDecays, &tbl, minBranching, excludeDecays);
      ImportDecays("Y", catalog, importedDecays, &tbl, minBranching, excludeDecays);
      
      // Instantiate new MDMSM object.
      auto pc = boost::make_shared<MDMSM>();
      
      // Populate annihilation channel list and add thresholds to threshold list.
      process_ann.resonances_thresholds.threshold_energy.push_back(2*mchi);
      auto channels = 
        daFunk::vec<string>("dbar_1, d_1", "ubar_2, u_2", "g, g", "dbar_3, d_3", "e+_1, e-_1", "ubar_1, u_1", "e+_2, e-_2", "Y, Y", "ubar_3, u_3", "dbar_2, d_2");
      auto p1 = 
        daFunk::vec<string>("dbar_1", "ubar_2", "g", "dbar_3", "e+_1", "ubar_1", "e+_2", "Y", "ubar_3", "dbar_2");
      auto p2 = 
        daFunk::vec<string>("d_1", "u_2", "g", "d_3", "e-_1", "u_1", "e-_2", "Y", "u_3", "d_2");
      
      for (unsigned int i = 0; i < channels.size(); ++i)
      {
        double mtot_final = 
        catalog.getParticleProperty(p1[i]).mass + 
        catalog.getParticleProperty(p2[i]).mass;  
        if (mchi*2 > mtot_final*0.5)
        {
          daFunk::Funk kinematicFunction = daFunk::funcM(pc, &MDMSM::sv, channels[i], tbl, 
          BEreq::CH_Sigma_V.pointer(), daFunk::var("v"));
          TH_Channel new_channel(daFunk::vec<string>(p1[i], p2[i]), kinematicFunction);
          process_ann.channelList.push_back(new_channel);
        }
        if (mchi*2 < mtot_final)
        {
          process_ann.resonances_thresholds.threshold_energy.
          push_back(mtot_final);
        }
      }
      
      if (spec.get(Par::Pole_Mass, "Y") >= 2*mchi) process_ann.resonances_thresholds.resonances.
          push_back(TH_Resonance(spec.get(Par::Pole_Mass, "Y"), tbl.at("Y").width_in_GeV));
      
      catalog.processList.push_back(process_ann);
      
      // Validate
      catalog.validate();
      
      result = catalog;
    } // function TH_ProcessCatalog
    
    void DarkMatter_ID_MDMSM(std::string& result){ result = "~chi"; }
    
    void DarkMatterConj_ID_MDMSM(std::string& result){ result = "~chi"; }
    
  } //namespace DarkBit
  
} //namespace Gambit

