//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  NREO model declaration. 
///
///  *********************************************
///
///  Authors (add name and date if you modify):
///   
///  \author Neal Avis Kozar
///  \date 2018 March
///
///  *********************************************


#ifndef __NREO_hpp__
#define __NREO_hpp__

#define MODEL NREO
	START_MODEL
	DEFINEPARS(j, m)
	DEFINEPARS(c0_1, c0_2, c0_3, c0_4, c0_5, c0_6, c0_7, c0_8, c0_9, c0_10, c0_11, c0_12, c0_13, c0_14, c0_15) 
	DEFINEPARS(c1_1, c1_2, c1_3, c1_4, c1_5, c1_6, c1_7, c1_8, c1_9, c1_10, c1_11, c1_12, c1_13, c1_14, c1_15)

	/// To declare this kind of relationship between a parameter 'my_par' and a 
	/// capability 'capability', one adds the following to the declaration of the 
	/// model containing 'my_par': MAP_TO_CAPABILITY(my_par, capability)
	MAP_TO_CAPABILITY(m, mwimpNREO)
	MAP_TO_CAPABILITY(j, jwimp)

	MAP_TO_CAPABILITY(c0_1, c0_1_cap)
	MAP_TO_CAPABILITY(c0_2, c0_2_cap)
	MAP_TO_CAPABILITY(c0_3, c0_3_cap)
	MAP_TO_CAPABILITY(c0_4, c0_4_cap)
	MAP_TO_CAPABILITY(c0_5, c0_5_cap)
	MAP_TO_CAPABILITY(c0_6, c0_6_cap)
	MAP_TO_CAPABILITY(c0_7, c0_7_cap)
	MAP_TO_CAPABILITY(c0_8, c0_8_cap)
	MAP_TO_CAPABILITY(c0_9, c0_9_cap)
	MAP_TO_CAPABILITY(c0_10, c0_10_cap)
	MAP_TO_CAPABILITY(c0_11, c0_11_cap)
	MAP_TO_CAPABILITY(c0_12, c0_12_cap)
	MAP_TO_CAPABILITY(c0_13, c0_13_cap)
	MAP_TO_CAPABILITY(c0_14, c0_14_cap)
	MAP_TO_CAPABILITY(c0_15, c0_15_cap)

	MAP_TO_CAPABILITY(c1_1, c1_1_cap)
	MAP_TO_CAPABILITY(c1_2, c1_2_cap)
	MAP_TO_CAPABILITY(c1_3, c1_3_cap)
	MAP_TO_CAPABILITY(c1_4, c1_4_cap)
	MAP_TO_CAPABILITY(c1_5, c1_5_cap)
	MAP_TO_CAPABILITY(c1_6, c1_6_cap)
	MAP_TO_CAPABILITY(c1_7, c1_7_cap)
	MAP_TO_CAPABILITY(c1_8, c1_8_cap)
	MAP_TO_CAPABILITY(c1_9, c1_9_cap)
	MAP_TO_CAPABILITY(c1_10, c1_10_cap)
	MAP_TO_CAPABILITY(c1_11, c1_11_cap)
	MAP_TO_CAPABILITY(c1_12, c1_12_cap)
	MAP_TO_CAPABILITY(c1_13, c1_13_cap)
	MAP_TO_CAPABILITY(c1_14, c1_14_cap)
	MAP_TO_CAPABILITY(c1_15, c1_15_cap)
#undef MODEL

#endif
