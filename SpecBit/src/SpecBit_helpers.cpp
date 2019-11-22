//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  Definitions of convenience (i.e. non-Gambit)
///  functions used by more than one SpecBit 
///  source file.
///
///  *********************************************
///
///  Authors (add name and date if you modify):
///
///  \author Tomas Gonzalo
///          (tomas.gonzalo@monash.edu)
///  \date 2019 Nov
///
///  *********************************************

#include "gambit/SpecBit/SpecBit_helpers.hpp"

namespace Gambit
{

  namespace SpecBit
  {

    /// @{

    /// Non-Gambit helper functions
    //  =======================================================================
    //  These are not known to Gambit, but perform helper tasks used by the
    //  Gambit module functions.

    /// Check that the SUSY spectrum has the canonical LSP for the model being scanned.
    void check_LSP(const Spectrum& spec, std::vector<int> LSPs)
    {
      // Fill the list of possible candidate LSPs is the are part of the spectrum
      // {~d, ~u, ~l, ~nu, ~g, ~chi0, ~chip}
      std::vector<int> LSP_candidates = {1000001,1000002,1000011,1000012,1000021,1000022,1000024};

      std::vector<double> m_LSP_candidates;
      for(int lsp_candidate : LSP_candidates)
        if(spec.has(Par::Pole_Mass, lsp_candidate, 0))
          m_LSP_candidates.push_back(std::abs(spec.get(Par::Pole_Mass, lsp_candidate, 0)));

      // Get the lightest of the canonical LSPs for this model
      double canonical_LSP = LSPs[0];
      for(int lsp : LSPs)
        if(spec.get(Par::Pole_Mass, canonical_LSP, 0) > spec.get(Par::Pole_Mass, lsp, 0))
          canonical_LSP = lsp;

      // Loop over candidates and invalidate if one is the lightest
      for(m_lsp_candidate : m_LSP_candidates)
      {
        if(spec.get(Par::Pole_Mass, canonical_LSP, 0) > m_lsp_candidate)
        {
          str canonical_LSP_name = Models::ParticleDB().long_name(canonical_LSP, 0);
          invalid_point().raise(canonical_LSP_name + " is not LSP.");
        }
      }

    }

    /// Helper to work with pointer
    void check_LSP(const Spectrum* spec, std::vector<int> LSPs)
    {
      check_LSP(*spec, LSPs);
    }

    /// Add gravitino mass to the spectrum and list of LSPs
    void add_gravitino_mass(Spectrum& spec, std::vector<int> &LSPs, double mG, const safe_ptr<Options>& runOptions)
    {
        static const double gmax = runOptions->getValueOrDef<double>(10.0, "max_permitted_gravitino_mass_GeV");
        if (mG > gmax)
        {
          std::ostringstream msg;
          msg << "Gravitino mass ("<<mG<<" GeV) is greater than permitted in *_lightgravitno models ("<<gmax<<" GeV).\n"
            << "If you know what you are doing(!), this behaviour can be modified with option max_permitted_gravitino_mass_GeV.";
          SpecBit_error().raise(LOCAL_INFO,msg.str());
        }

        spec.set(mG, Par::Pole_Mass, "~G");
        if (runOptions->getValueOrDef<bool>(true, "only_gravitino_LSP"))
          LSPs = {1000039};
        else
          LSPs.push_back(1000039);
    }


  }
}

