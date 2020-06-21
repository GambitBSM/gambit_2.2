//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  Models for QCD axions and axion-like particles.
///
///  *********************************************
///
///  Authors (add name and date if you modify):
///
///  \author Sebastian Hoof
///  \date 2016 Oct
///  \date 2017 Feb, May, Jul
///  \date 2018 Feb
///  \date 2019 Feb
///
///  *********************************************

#ifndef __GeneralALP_hpp__
#define __GeneralALP_hpp__

// General axion model with parametric temperature-dependent mass.
#define MODEL GeneralALP
  START_MODEL
  // Physical units: gagg [GeV^-1], gaee [dimensionless], gan [dimensionless]
  //                 gap [dimensionless], fa [GeV], ma0 [eV], Tchi [MeV],
  //                 beta [dimensionless], thetai [dimensionless]
  DEFINEPARS(gagg,gaee,gan,gap,fa,ma0,Tchi,beta,thetai)
#undef MODEL

// QCD axion model
#define MODEL QCDAxion
#define PARENT GeneralALP
  START_MODEL
  // Units for these parameters are the same as for the GeneralALP.
  DEFINEPARS(fa,Tchi,beta,thetai)
  // Physical units: LambdaChi [MeV], EoverN [dimensionless], CaggQCD [dimensionless]
  //                 Caee [dimensionless], Caee [dimensionless], Cap [dimensionless]
  DEFINEPARS(LambdaChi,EoverN,CaggQCD,Caee,Can,Cap)
  // Translation to parent, all defined in Axions.cpp:
  INTERPRET_AS_PARENT_FUNCTION(QCDAxion_to_GeneralALP)
#undef PARENT
#undef MODEL

// KSVZ axion model
#define MODEL KSVZAxion
#define PARENT QCDAxion
  START_MODEL
  // Units for these parameters are the same as for the QCDAxion.
  DEFINEPARS(fa,Tchi,beta,thetai,LambdaChi,EoverN,CaggQCD)
  // Physical units: Can0 [dimensionless], Cap0 [dimensionless], CaNTilde [dimensionless]
  DEFINEPARS(Can0,Cap0,CaNTilde)
  INTERPRET_AS_PARENT_FUNCTION(KSVZAxion_to_QCDAxion)
#undef PARENT
#undef MODEL

// DFSZ-I axion model
#define MODEL DFSZAxion_I
#define PARENT QCDAxion
  START_MODEL
  // Units for these parameters are the same as for the QCDAxion.
  DEFINEPARS(fa,Tchi,beta,thetai,LambdaChi,EoverN,CaggQCD)
  // Physical units: tanbeta [dimensionless], Can0 [dimensionless], Cap0 [dimensionless]
  //                 CaNTilde_uct [dimensionless], CaNTilde_dsb [dimensionless]
  DEFINEPARS(tanbeta,Can0,Cap0,CaNTilde_uct,CaNTilde_dsb)
  INTERPRET_AS_PARENT_FUNCTION(DFSZAxion_I_to_QCDAxion)
#undef PARENT
#undef MODEL

// DFSZ-II axion model
#define MODEL DFSZAxion_II
#define PARENT QCDAxion
  START_MODEL
  // Units for these parameters are the same as for the QCDAxion.
  DEFINEPARS(fa,Tchi,beta,thetai,LambdaChi,EoverN,CaggQCD)
  // Physical units: tanbeta [dimensionless], Can0 [dimensionless], Cap0 [dimensionless]
  //                 CaNTilde_uct [dimensionless], CaNTilde_dsb [dimensionless]
  DEFINEPARS(tanbeta,Can0,Cap0,CaNTilde_uct,CaNTilde_dsb)
  INTERPRET_AS_PARENT_FUNCTION(DFSZAxion_II_to_QCDAxion)
#undef PARENT
#undef MODEL

// ConstantMassALP model with temperature-independent mass and QCD-axion-inspired couplings
#define MODEL ConstantMassALP
#define PARENT GeneralALP
  START_MODEL
  // Units for these parameters are the same as for the GeneralALP.
  DEFINEPARS(fa,thetai)
  // Physical units: Cagg [dimensionless], Caee [dimensionless], Can [dimensionless], Cap [dimensionless], Lambda [MeV]
  DEFINEPARS(Cagg,Caee,Can,Cap,Lambda)
  INTERPRET_AS_PARENT_FUNCTION(ConstantMassALP_to_GeneralALP)
#undef PARENT
#undef MODEL

#endif
