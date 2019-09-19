//   GUM: GAMBIT Universal Models
//   **********************************
///  \file
///
///  Declarations of SARAH class
///
///  **********************************
///
///  \author Sanjay Bloor
///          (sanjay.bloor12@imperial.ac.uk)
///  \date 2017, 2018, 2019
///
///  \author Tomas Gonzalo
///          (tomas.gonzalo@monash.edu)
///  \date 2019 July, Aug
///
///  ***********************************

#ifndef SARAH_H
#define SARAH_H

#include "math_package.hpp"

namespace GUM
{
  class SARAH : public Math_Package
  {
    
    public:
  
      // SARAH constructor
      SARAH(std::string);

      // Load package
      void load_sarah();
    
      // Load model
      void load_model(std::string name);
    
      // Model checks
      bool check_model(std::string name);
    
      // Get model name
      std::string get_modelname();
    
      // Particle list
      void get_partlist(std::vector<Particle>&);
  
      // Parameters list
      void get_paramlist(std::vector<Parameter>&);

      // MINPAR, EXTPAR
      void get_minpar_extpar(std::vector<Parameter>&);

      // Get boundary conditions
      void get_boundary_conditions(std::vector<Parameter>&);

      // Add SPheno masses
      void add_SPheno_mass_names(std::vector<Particle>&);

      // Get tadpoles
      void get_tadpoles(std::vector<Parameter>&);
    
      // Outputs: CalcHEP, MadGraph, ...
      void write_ch_output();
      void write_micromegas_output();
      void write_madgraph_output();
      void write_spheno_output();
      void write_flexiblesusy_output();
      void write_vevacious_output();

  };

   
  // Everything
  void all_sarah(Options, std::vector<Particle>&, std::vector<Parameter>&, Outputs&, std::vector<std::string>&, std::map<std::string,bool>&);

} // namespace GUM
 
#endif
