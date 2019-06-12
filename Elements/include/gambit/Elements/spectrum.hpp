//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  Object for carrying around particle spectrum
///  data.
///  Primarily an interface to SLHAea objects,
///  carrying SLHA style information.
///
///  *********************************************
///
///  Authors:
///  <!-- add name and date if you modify -->
///
///  \author Ben Farmer
///          (benjamin.farmer@imperial.ac.uk)
///  \date 2019 June
//
///  *********************************************

#ifndef __Spectrum_hpp__
#define __Spectrum_hpp__

#include "gambit/Utils/util_functions.hpp"
#include "gambit/Utils/yaml_options.hpp"
#include "SLHAea/slhaea.h"

/// YAML overloads for mass cut and mass cut ratio constituents
namespace YAML
{

  typedef std::pair<std::string, std::pair<double, double> > sdd;
  typedef std::pair<std::pair<std::string,std::string>, std::pair<double, double> > ssdd;

  template<>
  struct convert<sdd>
  {
    static Node encode(const sdd& rhs)
    {
      Node node;
      node.push_back(rhs.first);
      node.push_back(rhs.second.first);
      node.push_back(rhs.second.second);
      return node;
    }

    static bool decode(const Node& node, sdd& rhs)
    {
      if(!node.IsSequence() || node.size() != 3) return false;
      rhs.first         = node[0].as<std::string>();
      rhs.second.first  = node[1].as<double>();
      rhs.second.second = node[2].as<double>();
      return true;
    }
  };

  template<>
  struct convert<ssdd>
  {
    static Node encode(const ssdd& rhs)
    {
      Node node;
      node.push_back(rhs.first.first);
      node.push_back(rhs.first.second);
      node.push_back(rhs.second.first);
      node.push_back(rhs.second.second);
      return node;
    }

    static bool decode(const Node& node, ssdd& rhs)
    {
      if(!node.IsSequence() || node.size() != 4) return false;
      rhs.first.first   = node[0].as<std::string>();
      rhs.first.second  = node[1].as<std::string>();
      rhs.second.first  = node[2].as<double>();
      rhs.second.second = node[3].as<double>();
      return true;
    }
  };

}


namespace Gambit
{

   /// Less confusing name for SLHAea container class
   typedef SLHAea::Coll SLHAstruct;

   /// Class for interfacing to spectrum information
   class Spectrum
   {
      public:

         /// Typedefs for making it easier to manipulate mass cut and mass ratio cut info.
         /// @{
         typedef std::vector<YAML::sdd>  mc_info;
         typedef std::vector<YAML::ssdd> mr_info;
         /// @}

      private:

         /// Variables
         /// @{
         mc_info mass_cuts;
         mr_info mass_ratio_cuts;
         /// @}

         ///Calculate Wolfenstein rho+i*eta from rhobar and etabar
         static std::complex<double> rhoplusieta(double, double, double, double);

         /// Wrapped SLHAea object
         SLHAstruct mySLHAea;

         /// Contents requirements for this spectrum
         SpectrumContents::Contents myContents;

      public:

         /// @{ Constructors/Destructors
         /// Need custom copy and move constructors plus copy-assignment operator
         /// in order to manage the unique_ptrs properly.

         /// Default constructor
         Spectrum();

         /// Construct from SLHAea object (also specifying what SpectrumContents should apply, which defines how to interpret the SLHAea blocks)
         Spectrum(const SLHAstruct& slha, const SpectrumContents::Contents& contents);

         /// Set constraints on masses and mass ratios that cause the spectrum to be declared "invalid" if they are violated
         void set_mass_cuts(const mc_info&);
         void set_mass_ratio_cuts(const mr_info&);
         /// @}

         /// Check the that the spectrum satisifies any mass cuts requested from the yaml file.
         void check_mass_cuts();

         /// @{ Pole mass getters
         //
         //
         // bool auto_check_antiparticle_name
         
         /// Master checker function to see if a parameter request matches the spectrum contents
         bool has(const Par::Tags partype, std::string& name, std::vector<int> indices, bool auto_check_antiparticle_name=true) const;
         
         /// Getters to access pole masses in spectrum 
         bool   has(const Par::Tags partype, const std::string& mass) const;
         double get(const Par::Tags partype, const std::string& mass) const;
         bool   has(const Par::Tags partype, const std::string& mass, const int index) const;
         double get(const Par::Tags partype, const std::string& mass, const int index) const;

         /// @{ PDB getter/checker overloads
         bool   has(const Par::Tags partype, const int pdg_code, const int context) const;
         double get(const Par::Tags partype, const int pdg_code, const int context) const;
         bool   has(const Par::Tags partype, const std::pair<int,int> pdgpr) const;
         double get(const Par::Tags partype, const std::pair<int,int> pdgpr) const;
         bool   has(const Par::Tags partype, const std::pair<str,int> shortpr) const;
         double get(const Par::Tags partype, const std::pair<str,int> shortpr) const;
         /// @}

         /// @{ Getters which first check the sanity of the thing they are returning
         double safeget(const Par::Tags partype, const std::string& mass) const;
         double safeget(const Par::Tags partype, const std::string& mass, const int index) const;
         double safeget(const Par::Tags partype, const int pdg_code, const int context) const;
         double safeget(const Par::Tags partype, const std::pair<int,int> pdgpr) const;
         double safeget(const Par::Tags partype, const std::pair<str,int> shortpr) const;
         /// @}

         /* Setter declarations, for manually overwriting parameter values (directly changes wrapped SLHAea object contents)
            Note; these are NON-CONST */
         void set(const Par::Tags, const double, const str&);
         void set(const Par::Tags, const double, const str&, const int);
         void set(const Par::Tags, const double, const str&, const int, const int);

         /* Setters for setting values of many parameters at once, by iterating over the supplied string names or indices, or both */
         void set_many(const Par::Tags, const double, const std::vector<str>&);
         void set_many(const Par::Tags, const double, const std::vector<str>&, const std::vector<int>);
         void set_many(const Par::Tags, const double, const std::vector<str>&, const int);
         void set_many(const Par::Tags, const double, const str&, const std::vector<int>);

         /// @}

         /// SLHAea object getter
         /// Return underlying SLHAea object (with conversions between SLHA1/2 definitions if applicable)
         SLHAstruct getSLHAea(int) const;

         /// Output spectrum contents as an SLHA file, using getSLHAea.
         void writeSLHAfile(int, const str&) const;

         /// Helper function to drop SLHA files
         void drop_SLHAs_if_requested(const safe_ptr<Options>&, const str&);

         /// PDG code translation map, for special cases where an SLHA file has been read in and the PDG codes changed.
         const std::map<int, int>& PDG_translator() const;

         /// CKM Wolfenstein (lambda, A, rhobar, etabar) --> V_qq standard parameterisation convertors
         /// @{
         static double Wolf2V_ud(double, double, double, double);
         static double Wolf2V_us(double, double, double, double);
         static std::complex<double> Wolf2V_ub(double, double, double, double);
         static std::complex<double> Wolf2V_cd(double, double, double, double);
         static std::complex<double> Wolf2V_cs(double, double, double, double);
         static double Wolf2V_cb(double, double, double, double);
         static std::complex<double> Wolf2V_td(double, double, double, double);
         static std::complex<double> Wolf2V_ts(double, double, double, double);
         static double Wolf2V_tb(double, double, double, double);
         /// @}
   };

} // end namespace Gambit


#endif
