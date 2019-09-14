//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  Functions for analyses that use interpolated yields.
///
///  *********************************************
///
///  Authors (add name and date if you modify):
///
///  \author Martin White
///          (martin.white@adelaide.edu.au)
///  \author Andre Scaffidi
///          (andre.scaffidi@adelaide.edu.au)
///  \date 2019 Aug
///
///  Notes:
///
///   - has been put together for the DMEFT project
///   - could probably introduce a better capability
///     structure if we decide to use the functionality
///     for other models
///
///  *********************************************

#include <cmath>
#include <string>
#include <iostream>
#include <fstream>
#include <memory>
#include <numeric>
#include <sstream>
#include <vector> 
#include <iomanip>
#include <math.h>
#include <cstring>
#include <stdlib.h>
using namespace std;



#include "gambit/ColliderBit/analyses/Analysis.hpp"
#include "gambit/Elements/gambit_module_headers.hpp"
#include "gambit/ColliderBit/ColliderBit_rollcall.hpp"
#include "gambit/Utils/ascii_table_reader.hpp"

// Needs GSL 2 
#include <gsl/gsl_math.h>
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>
#include <gsl/gsl_sf_gamma.h>

#include "Eigen/Eigenvalues"
#include "Eigen/Eigen"


namespace Gambit
{

  namespace ColliderBit
  {
    
// ---------------------------------- FUNCTION OF INTEREST ------------------------------------------------

    // ---- Define my interpolation fucntions here as well as get data functions
  const char* colliderbitdata_path = GAMBIT_DIR "/ColliderBit/data/"; 
  #define PI 3.14159265

  // Initialize all data
  static const size_t data_INC           = 20;
  static const size_t data_SIZE          = pow(data_INC,2);
  static const size_t cms_bin_size       = 22;
  static const size_t atlas_bin_size     = 11;
  double THETA[data_SIZE];
  double MASS[data_SIZE];
  double nJets[data_SIZE];
  double CS[data_SIZE];
  double Total_events[data_SIZE];
  double MET_HIST_CMS[data_SIZE][cms_bin_size];
  double MET_HIST_ATLAS[data_SIZE][atlas_bin_size];

  // Define just mass and angle arrays
  double theta[data_INC];
  double mass[data_INC];

  double BilinearInterpolation(double q11, double q12, double q21, double q22, 
      double x1, double x2, double y1, double y2, double x, double y) 
  {
    double x2x1, y2y1, x2x, y2y, yy1, xx1;
    x2x1 = x2 - x1;
    y2y1 = y2 - y1;
    x2x = x2 - x;
    y2y = y2 - y;
    yy1 = y - y1;
    xx1 = x - x1;
    return 1.0 / (x2x1 * y2y1) * (
      q11 * x2x * y2y +
      q21 * xx1 * y2y +
      q12 * x2x * yy1 +
      q22 * xx1 * yy1
    );
  }



  void Get_data(const char* file_grid,const char* file_X_Y, const char* file_met_hist ){
    // Load Data

    ifstream mb(file_met_hist);
    // cout << file_grid << " " << file_X_Y << " " << file_met_hist <<endl;

    if (file_met_hist == colliderbitdata_path+string("DMEFT/met_hist_ATLAS_C62_C63.txt") || file_met_hist ==colliderbitdata_path+string("DMEFT/met_hist_ATLAS_C61_C64.txt")){

      // std::string str;

      // while(getline(mb,str)){

      for(int row = 0; row < data_SIZE; row++) {  
        for(int column = 0; column < atlas_bin_size; column++){
          mb >> MET_HIST_ATLAS[row][column];
        }
        }
      // }

    }		
    if (file_met_hist == colliderbitdata_path+string("DMEFT/met_hist_CMS_C62_C63.txt") || file_met_hist ==colliderbitdata_path+string("DMEFT/met_hist_CMS_C61_C64.txt")){
      // std::string str;
      // while(getline(mb,str)){
        for(int row = 0; row < data_SIZE; row++) {  
          for(int column = 0; column < cms_bin_size; column++){
            mb >> MET_HIST_CMS[row][column];
            // cout <<  MET_HIST_CMS[row][column]<< " "<< row << " "<< column<<endl;
          }
        }
      // }
    }	

    // Rest of grid

    float p1,p2,p3,p4,p5;
    FILE * pp = fopen(file_grid,"r");              // file containing numbers in 5 columns 
    for (int ll = 0; ll < data_SIZE; ++ll){
      fscanf(pp,"%f %f %f %f", &p1,&p2,&p3,&p4);
      // cout << p1 << endl;
      MASS[ll] = p1; 
      // cout <<"Mass_check = "<< MASS[ll]<<" "<<ll<<endl;
      THETA[ll]= p2;
      nJets[ll]= p3;
      CS[ll]   = p4;   
    }
    fclose(pp);


    float var1,var2;
    FILE * fp = fopen(file_X_Y,"r");             
    for (int ll = 0; ll < data_INC; ++ll){
      fscanf(fp,"%f %f", &var1, &var2);
      mass[ll] = var1;
      theta[ll]= var2;
      }
    fclose(fp);	


  }
  // ---------------------------------------------------- //
  // Interpolation functions // 
  // ---------------------------------------------------- //

  double * Acceptance_CS(float m,float O1,float O2, const char* pair, const char* experiment){
    char const *i = (colliderbitdata_path+string("DMEFT/met_hist_ATLAS_C61_C64.txt")).c_str();
    char const *j = (colliderbitdata_path+string("DMEFT/met_hist_ATLAS_C62_C63.txt")).c_str();
    char const *k = (colliderbitdata_path+string("DMEFT/met_hist_CMS_C61_C64.txt")).c_str();
    char const *l = (colliderbitdata_path+string("DMEFT/met_hist_CMS_C62_C63.txt")).c_str();
    char const *b = (colliderbitdata_path+string("DMEFT/X_Y_ATLAS_C62_C63.txt")).c_str();
    char const *a = (colliderbitdata_path+string("DMEFT/grid_output_ATLAS_C62_C63.txt")).c_str();
    char const *h = (colliderbitdata_path+string("DMEFT/X_Y_CMS_C61_C64.txt")).c_str();
    char const *g = (colliderbitdata_path+string("DMEFT/grid_output_CMS_C61_C64.txt")).c_str();
    char const *f = (colliderbitdata_path+string("DMEFT/X_Y_ATLAS_C61_C64.txt")).c_str();
    char const *e = (colliderbitdata_path+string("DMEFT/grid_output_ATLAS_C61_C64.txt")).c_str();
    char const *d = (colliderbitdata_path+string("DMEFT/X_Y_CMS_C62_C63.txt")).c_str();
    char const *c = (colliderbitdata_path+string("DMEFT/grid_output_CMS_C62_C63.txt")).c_str();
    int met_bin_size;
    double ** MET_HIST = new double*[data_SIZE];

    
    if (experiment=="ATLAS" && pair == "23"){
      Get_data(a,b,j);
      met_bin_size = atlas_bin_size;

      // double** MET_HIST = new double*[data_SIZE];
      for(int i = 0; i < data_SIZE; ++i)
        MET_HIST[i] = new double[met_bin_size];

      // Assign met histogram to current experiment
      for (int kk = 0; kk<data_SIZE;++kk){
        for (int j = 0; j<met_bin_size;++j){
          MET_HIST[kk][j] = MET_HIST_ATLAS[kk][j];
          }
        }
    }

    else if (experiment=="CMS" && pair == "23"){
      Get_data(c,d,l);
      met_bin_size = cms_bin_size;
      
      for(int i = 0; i <= data_SIZE; ++i){
        MET_HIST[i] = new double[met_bin_size];
      }
      // Assign met histogram to current experiment
      for (int kk = 0; kk<data_SIZE;++kk){
        for (int j = 0; j<met_bin_size;++j){
          MET_HIST[kk][j] = MET_HIST_CMS[kk][j];
          }
        }
      
    }

    else if (experiment=="ATLAS" && pair == "14"){
      Get_data(e,f,i);
      met_bin_size = atlas_bin_size;
      // double** MET_HIST = new double*[data_SIZE];
      for(int i = 0; i < data_SIZE; ++i)
        MET_HIST[i] = new double[met_bin_size];
      // Assign met histogram to current experiment
      for (int kk = 0; kk<data_SIZE;++kk){
        for (int j = 0; j<met_bin_size;++j){
          MET_HIST[kk][j] = MET_HIST_ATLAS[kk][j];
          }
        }
    }

    else if (experiment=="CMS" && pair == "14"){
      Get_data(g,h,k);
      met_bin_size = cms_bin_size;
      // double** MET_HIST = new double*[data_SIZE];
      for(int i = 0; i < data_SIZE; ++i)
        MET_HIST[i] = new double[met_bin_size];
      // Assign met histogram to current experiment
      for (int kk = 0; kk<data_SIZE;++kk){
        for (int j = 0; j<met_bin_size;++j){
          MET_HIST[kk][j] = MET_HIST_CMS[kk][j];
          }
        }
    }


    // DEBUG stuff
    // cout << "Met bin size" << " "<< met_bin_size<<" "<<data_SIZE<<endl;


    // for (int kk = 0; kk<data_SIZE;++kk){
    // 	for (int j = 0; j<met_bin_size-1;++j){
    // 		// cout << kk << " "<< j<<" "<< MET_HIST[kk][j] <<endl;
    // 		if (MET_HIST[kk][j] == 0){
    // 			cout << "Zero! Bin = " <<" "<<j<<": mass,theta ="<<MASS[kk]<< " "<<THETA[kk]<< " " << "CMS val = "<< " "<<MET_HIST_CMS[kk][j] << endl;
    // 		}
    // 	}
    // }
    // cout << "Done!"<<endl;



    
    // Get scale factor and phase theta
    
    double Norm,th;

    if (O1==0){
      Norm = pow(O2,2);
      // cout << " O1 is zero"<< endl;
    }
    else if (O2==0){
      Norm = pow(O1,2);
      // cout << " O2 is zero"<< endl;
    }
    else{
      th    = 0.5*asin(2*O1*O2/(pow(O1,2)+pow(O2,2)));
      if (O1*O2 < 0){
        th = th + PI;
      }
      cout << "Theta = "<< th<<endl;
      Norm  = 2*O1*O2/(sin(2.0*th));
    }

    // Checks to go ahead with interpolation

    if (m<mass[0] || m>mass[data_INC-1]){
      cout<<" Error! Mass param out of range with value "  << m <<  " Exiting..."<<endl;
      std::exit(EXIT_SUCCESS);
      }
    else if (th<theta[0] || th>theta[data_INC-1]){
      cout<<" Error! Theta param out of range with value " << th << " Exiting..."<<endl;
      std::exit(EXIT_SUCCESS);
      }

    cout << "MJW a" << endl;
    
    // Get x1,2 y1,2 : Mass and theta coordinates for interpolation
    double x1,x2,y1,y2;
    for(int ii = 0; ii < data_INC-1; ++ii) {
      if (m >= mass[ii] && m <= mass[ii+1]){
        x1 = mass[ii];
        x2 = mass[ii+1];
        break;
        }
      }
    for(int jj = 0; jj < data_INC-1; ++jj) {
      if (th >= theta[jj] && th <= theta[jj+1]){
        y1 = theta[jj];
        y2 = theta[jj+1];
        break;
        }
      }

    // Get C's
    cout << "MJW b" << endl;
    double C11=0.0 ,C12=0.0,C21=0.0,C22=0.0;

    // Define Q's as array: One Q type for each met bin.

    double Q11[met_bin_size]={0} ,Q12[met_bin_size]={0},Q21[met_bin_size]={0},Q22[met_bin_size]={0};

    double* accep = new double[met_bin_size]; 
    
    cout << "MJW c" << endl;
    // NJets and Cross-section
    for (int Emiss = 0; Emiss < met_bin_size-1; ++Emiss ) {
      while (Q11[Emiss]==0.0 || Q12[Emiss]==0.0 || Q21[Emiss]== 0.0 || Q22[Emiss]==0.0 || C11==0.0 || C12==0.0 || C21== 0.0 || C22==0.0){ 
        for(int kk = 0; kk < data_SIZE; ++kk) {
          // cout << MASS[kk]<<" "<< THETA[kk]<< "|  |"<<MET_HIST[kk][Emiss]<<" " << kk<< " |     |" << x2<<" " << y2 <<" "<< Q22[Emiss]<<endl;
          
          if (MASS[kk]==x1 && THETA[kk]==y1){
            // Q11[Emiss] = nJets[kk];
            Q11[Emiss] = MET_HIST[kk][Emiss];
            C11 = CS[kk];

            }
          else if (MASS[kk]==x1 && THETA[kk]==y2){
            Q12[Emiss] = MET_HIST[kk][Emiss];
            C12 = CS[kk];

            }
          else if (MASS[kk]==x2 && THETA[kk]==y1){
            Q21[Emiss] = MET_HIST[kk][Emiss];
            C21 = CS[kk];

            }
          else if (MASS[kk]==x2 && THETA[kk]==y2){
            Q22[Emiss] = MET_HIST[kk][Emiss];
            C22 = CS[kk];

            // cout << " Q22"<< " "<< Q22[Emiss]<<endl;

            }	
          }	
        }

      cout << "MJW d" << endl;

      // Need to check if this is good with Martin.
      double res =  Norm*BilinearInterpolation(Q11[Emiss], Q12[Emiss], Q21[Emiss], Q22[Emiss], x1, x2, y1, y2, m, th)*Norm*BilinearInterpolation(C11, C12, C21, C22, x1, x2, y1, y2, m, th); 
    

      accep[Emiss] = res;
    }

    cout << "Returning" << endl;
    
    return accep;

    cout << "Still here" << endl;
    
    std::fill_n(THETA,data_SIZE,0);
    std::fill_n(MASS,data_SIZE,0);
    std::fill_n(nJets,data_SIZE,0);
    std::fill_n(CS,data_SIZE,0);
    std::fill_n(Total_events,data_SIZE,0);
    MET_HIST_CMS[data_SIZE][cms_bin_size] = {};
    MET_HIST_ATLAS[data_SIZE][atlas_bin_size] = {};
    MET_HIST[data_SIZE][met_bin_size] = {};

  }

double *  Acc_Eff_CS(float m,float C61,float C62,float C63, float C64 , const char* exper_){
	char const *tt = "23";
	char const *of = "14";
	int met_bin_size;

	if (exper_=="ATLAS"){
		met_bin_size = atlas_bin_size;
	}
	else if (exper_=="CMS"){
		met_bin_size = cms_bin_size;
	}

	double* YIELDS = new double[met_bin_size]; 
	
	double* A23;
	double* A14;

	A23 = Acceptance_CS(m,C62,C63,tt,exper_);
	A14 = Acceptance_CS(m,C61,C64,of,exper_);

	for (int ii = 0; ii < met_bin_size-1; ++ii){
		YIELDS[ii] = A23[ii] + A14[ii];
		// YIELDS[ii] = A14[ii];

	}

	return YIELDS;
}

    void DMEFT_results(AnalysisDataPointers &result){  


      // This routine will get the yields for both the ATLAS and CMS monojet analyses
      // The results are stored in a vector of AnalysisData objects, which includes backgrounds yields, uncertainties and correlations
      
      // Andre: have put a dummy file of cross-section data here
      // Will need to be replaced by the relevant file

      // Also you will need to add contributions from the other intefering operators

      // Am assuming for now that this is fed the actual Wilson coefficients
      // You will need to add code that maps these to the mixing angle, etc
      
    
      float C61 = *Pipes::DMEFT_results::Param["C61"];
      float C62 = *Pipes::DMEFT_results::Param["C62"];
      float C63 = *Pipes::DMEFT_results::Param["C63"];
      float C64 = *Pipes::DMEFT_results::Param["C64"];
            
      // ***** What about DM mass?


      // Andre: will need too add interpolators for each bin (or some smarter way to do it for all bins and store the results)
      // ColliderBitInterpolator2D cross_C61_C64(colliderbitdata_path+"DMEFT/test_crosssec.dat","bicubic");
      // ColliderBitInterpolator2D eff_C61_C64_ATLAS(colliderbitdata_path+"DMEFT/test_eff.dat","bicubic");
      //ColliderBitInterpolator2D eff_C61_C64_CMS(...)

      // double *yield_C61_C64_ATLAS;

      // Andre: for now am assuming that the yield is just 1000. * eff * crosssec
      // This will need to be updated
      // // Check that the interpolators are valid
      // if((C61 < cross_C61_C64.lower_x()) ||
      //   (C61 > cross_C61_C64.upper_x()) ||
      //   (C64 < cross_C61_C64.lower_y()) ||
      //   (C64 > cross_C61_C64.upper_y())){
        
	    //   td::cerr << "Interpolator is out of bounds in ColliderBit DMEFT calculations. The results will be totally wrong." << std::endl;
	
      //    }

      // else {
      //   //Andre to replace this calculation
      //   yield_C61_C64_ATLAS = 1000. * cross_C61_C64.interpolate(C61,C64) * eff_C61_C64_ATLAS.interpolate(C61,C64);
      //   //yield_C61_C64_CMS = 1000. * cross_C61_C64.interpolate(C61,C64) * eff_C61_C64_CMS.interpolate(C61,C64);
      // }
	  
      // std::cout << "Yield: " << yield_C61_C64_ATLAS << std::endl;

      // Put CMS signal region data into an AnalysisData object
      // Note: includes covariance matrix
      // static const size_t NUMSR = 22;

      // Will need to set the vector _srnums to hold the interpolated yields in each bin

      // Andre needs to put signal numbers for CMS bins here (output from interpolator)


        // ----------------------------------------------------------------------------------------------------//
         // --------------------------------CMS---------------------------------------------------------//

      // Test the function to see if it compiles. 
      double mass = 150;
	    
      double *_srnums;

      _srnums = Acc_Eff_CS(mass,C61,C62,C63,C64,"CMS");
      

      static const double OBSNUM[cms_bin_size] = {
                              136865, 74340, 42540, 25316, 15653, 10092, 8298, 4906, 2987, 2032, 1514,
                              926, 557, 316, 233, 172, 101, 65, 46, 26, 31, 29};
      static const double BKGNUM[cms_bin_size] = {
                                          134500, 73400, 42320, 25490, 15430, 10160, 8480, 4865, 2970, 1915, 1506,
                                          844, 526, 325, 223, 169, 107, 88.1, 52.8, 25.0, 25.5, 26.9
                                             };
      static const double BKGERR[cms_bin_size] = { 
                                          	3700, 2000, 810, 490, 310, 170, 140, 95, 49, 33, 32, 18, 14, 12, 9, 8, 6, 5.3, 3.9, 2.5, 2.6, 2.8
                                                };

      std::vector<SignalRegionData> cmsBinnedResults;
      
      for (size_t ibin = 0; ibin < cms_bin_size; ++ibin) {
        	std::stringstream ss; ss << "sr-" << ibin;
	        cmsBinnedResults.push_back(SignalRegionData(ss.str(), OBSNUM[ibin], {_srnums[ibin],  0.}, {BKGNUM[ibin], BKGERR[ibin]}));
      }
      
      static const std::vector< std::vector<double> > BKGCOV = {
	{  1.37e+07,  7.18e+06,  2.58e+06,  1.54e+06,  9.29e+05,  4.28e+05,  3.26e+05,  2.04e+05,  8.34e+04,  5.37e+04,  4.62e+04,  2.33e+04,  1.45e+04,  1.20e+04,  6.66e+03,  7.99e+03,  4.00e+03,  1.57e+03,  0.00e+00,  1.30e+03,  3.85e+02, -4.14e+02 },
	{  7.18e+06,  4.00e+06,  1.38e+06,  8.43e+05,  5.02e+05,  2.28e+05,  1.74e+05,  1.05e+05,  4.51e+04,  2.84e+04,  2.30e+04,  1.22e+04,  7.56e+03,  6.48e+03,  3.24e+03,  4.00e+03,  2.28e+03,  1.06e+03,  1.56e+02,  8.00e+02,  3.64e+02, -1.68e+02 },
	{  2.58e+06,  1.38e+06,  6.56e+05,  3.57e+05,  2.18e+05,  1.07e+05,  8.73e+04,  5.31e+04,  2.34e+04,  1.50e+04,  1.35e+04,  7.00e+03,  4.20e+03,  3.30e+03,  2.26e+03,  1.81e+03,  1.12e+03,  6.44e+02,  2.21e+02,  3.04e+02,  1.47e+02,  2.27e+01 },
	{  1.54e+06,  8.43e+05,  3.57e+05,  2.40e+05,  1.32e+05,  6.58e+04,  5.14e+04,  3.17e+04,  1.44e+04,  9.22e+03,  8.15e+03,  4.06e+03,  2.88e+03,  2.00e+03,  1.32e+03,  1.25e+03,  7.06e+02,  3.64e+02,  5.73e+01,  1.59e+02,  7.64e+01, -2.74e+01 },
	{  9.29e+05,  5.02e+05,  2.18e+05,  1.32e+05,  9.61e+04,  4.11e+04,  3.21e+04,  1.88e+04,  8.81e+03,  5.73e+03,  5.46e+03,  2.57e+03,  1.78e+03,  1.34e+03,  6.98e+02,  9.18e+02,  4.28e+02,  1.64e+02,  3.63e+01,  1.32e+02,  1.05e+02, -8.68e+00 },
	{  4.28e+05,  2.28e+05,  1.07e+05,  6.58e+04,  4.11e+04,  2.89e+04,  1.76e+04,  1.07e+04,  5.16e+03,  2.92e+03,  2.83e+03,  1.62e+03,  9.76e+02,  8.77e+02,  3.82e+02,  4.49e+02,  2.04e+02,  1.08e+02,  9.94e+01,  1.02e+02,  3.98e+01,  4.76e+00 },
	{  3.26e+05,  1.74e+05,  8.73e+04,  5.14e+04,  3.21e+04,  1.76e+04,  1.96e+04,  9.18e+03,  4.39e+03,  2.82e+03,  2.46e+03,  1.39e+03,  9.21e+02,  7.39e+02,  5.17e+02,  3.70e+02,  2.35e+02,  9.65e+01,  8.19e+01,  4.20e+01,  1.82e+01,  3.14e+01 },
	{  2.04e+05,  1.04e+05,  5.31e+04,  3.17e+04,  1.88e+04,  1.07e+04,  9.18e+03,  9.02e+03,  2.61e+03,  1.72e+03,  1.70e+03,  8.55e+02,  4.52e+02,  4.67e+02,  2.48e+02,  2.66e+02,  1.54e+02,  5.04e+01,  3.33e+01,  1.19e+01,  3.21e+01,  7.98e+00 },
	{  8.34e+04,  4.51e+04,  2.34e+04,  1.44e+04,  8.81e+03,  5.16e+03,  4.39e+03,  2.61e+03,  2.40e+03,  9.22e+02,  8.94e+02,  4.67e+02,  2.13e+02,  2.41e+02,  1.41e+02,  1.29e+02,  4.70e+01,  4.41e+01,  7.64e+00,  2.08e+01,  2.55e+01,  5.49e+00 },
	{  5.37e+04,  2.84e+04,  1.50e+04,  9.22e+03,  5.73e+03,  2.92e+03,  2.82e+03,  1.72e+03,  9.22e+02,  1.09e+03,  5.17e+02,  3.03e+02,  1.62e+02,  1.47e+02,  8.91e+01,  8.18e+01,  3.17e+01,  2.10e+01,  1.29e+00,  7.42e+00,  7.72e+00,  4.62e+00 },
	{  4.62e+04,  2.30e+04,  1.35e+04,  8.15e+03,  5.46e+03,  2.83e+03,  2.46e+03,  1.70e+03,  8.94e+02,  5.17e+02,  1.02e+03,  2.65e+02,  1.57e+02,  1.61e+02,  9.22e+01,  7.94e+01,  3.84e+01,  3.39e+00, -1.25e+00,  1.44e+01,  3.33e+00, -8.96e-01 },
	{  2.33e+04,  1.22e+04,  7.00e+03,  4.06e+03,  2.57e+03,  1.62e+03,  1.39e+03,  8.55e+02,  4.67e+02,  3.03e+02,  2.65e+02,  3.24e+02,  8.57e+01,  9.07e+01,  5.83e+01,  3.02e+01,  2.70e+01,  2.00e+01,  7.02e+00,  2.25e+00,  5.15e+00,  7.06e+00 },
	{  1.45e+04,  7.56e+03,  4.20e+03,  2.88e+03,  1.78e+03,  9.76e+02,  9.21e+02,  4.52e+02,  2.13e+02,  1.62e+02,  1.57e+02,  8.57e+01,  1.96e+02,  5.21e+01,  3.91e+01,  3.92e+01,  2.69e+01,  8.90e+00,  6.55e+00,  0.00e+00,  1.46e+00,  1.57e+00 },
	{  1.20e+04,  6.48e+03,  3.30e+03,  2.00e+03,  1.34e+03,  8.77e+02,  7.39e+02,  4.67e+02,  2.41e+02,  1.47e+02,  1.61e+02,  9.07e+01,  5.21e+01,  1.44e+02,  3.02e+01,  2.02e+01,  1.44e+01,  3.18e+00,  4.68e-01,  4.50e+00,  2.18e+00,  3.02e+00 },
	{  6.66e+03,  3.24e+03,  2.26e+03,  1.32e+03,  6.98e+02,  3.82e+02,  5.17e+02,  2.48e+02,  1.41e+02,  8.91e+01,  9.22e+01,  5.83e+01,  3.91e+01,  3.02e+01,  8.10e+01,  1.15e+01,  1.19e+01,  7.63e+00,  3.16e+00, -2.25e-01,  1.40e+00,  2.52e+00 },
	{  7.99e+03,  4.00e+03,  1.81e+03,  1.25e+03,  9.18e+02,  4.49e+02,  3.70e+02,  2.66e+02,  1.29e+02,  8.18e+01,  7.94e+01,  3.02e+01,  3.92e+01,  2.02e+01,  1.15e+01,  6.40e+01,  1.92e+00, -1.27e+00, -3.12e-01,  1.40e+00,  2.70e+00, -6.72e-01 },
	{  4.00e+03,  2.28e+03,  1.12e+03,  7.06e+02,  4.28e+02,  2.04e+02,  2.35e+02,  1.54e+02,  4.70e+01,  3.17e+01,  3.84e+01,  2.70e+01,  2.69e+01,  1.44e+01,  1.19e+01,  1.92e+00,  3.60e+01,  5.09e+00,  3.74e+00, -1.65e+00,  1.40e+00,  1.51e+00 },
	{  1.57e+03,  1.06e+03,  6.44e+02,  3.64e+02,  1.64e+02,  1.08e+02,  9.65e+01,  5.04e+01,  4.41e+01,  2.10e+01,  3.39e+00,  2.00e+01,  8.90e+00,  3.18e+00,  7.63e+00, -1.27e+00,  5.09e+00,  2.81e+01,  6.20e-01, -1.19e+00,  5.51e-01, -4.45e-01 },
	{  0.00e+00,  1.56e+02,  2.21e+02,  5.73e+01,  3.63e+01,  9.95e+01,  8.19e+01,  3.33e+01,  7.64e+00,  1.29e+00, -1.25e+00,  7.02e+00,  6.55e+00,  4.68e-01,  3.16e+00, -3.12e-01,  3.74e+00,  6.20e-01,  1.52e+01,  7.80e-01,  3.04e-01,  1.64e+00 },
	{  1.30e+03,  8.00e+02,  3.04e+02,  1.59e+02,  1.32e+02,  1.02e+02,  4.20e+01,  1.19e+01,  2.08e+01,  7.42e+00,  1.44e+01,  2.25e+00,  0.00e+00,  4.50e+00, -2.25e-01,  1.40e+00, -1.65e+00, -1.19e+00,  7.80e-01,  6.25e+00,  1.30e-01,  6.30e-01 },
	{  3.85e+02,  3.64e+02,  1.47e+02,  7.64e+01,  1.05e+02,  3.98e+01,  1.82e+01,  3.21e+01,  2.55e+01,  7.72e+00,  3.33e+00,  5.15e+00,  1.46e+00,  2.18e+00,  1.40e+00,  2.70e+00,  1.40e+00,  5.51e-01,  3.04e-01,  1.30e-01,  6.76e+00,  5.82e-01 },
	{ -4.14e+02, -1.68e+02,  2.27e+01, -2.74e+01, -8.68e+00,  4.76e+00,  3.14e+01,  7.98e+00,  5.49e+00,  4.62e+00, -8.96e-01,  7.06e+00,  1.57e+00,  3.02e+00,  2.52e+00, -6.72e-01,  1.51e+00, -4.45e-01,  1.64e+00,  6.30e-01,  5.82e-01,  7.84e+00 }
      };

      Eigen::MatrixXd m_BKGCOV(22,22);
      for (int i = 0; i < 22; i++)
	    m_BKGCOV.row(i) = Eigen::VectorXd::Map(&BKGCOV[i][0],BKGCOV[i].size());
      
      AnalysisData  * cmsData = new AnalysisData(cmsBinnedResults, m_BKGCOV);
      cmsData->analysis_name = "CMS_13TeV_MONOJET_36invfb_interpolated";

  // ----------------------------------------------------------------------------------------------------//
  // ---------------------------------ATLAS----------------------------------------------------------//

      // Now put the ATLAS data into an equivalent object
      // Andre to add the relevant lines
      

      std::cout << "Making signal numbers" << std::endl;
      double *_srnums_ATLAS;

      _srnums_ATLAS = Acc_Eff_CS(mass,C61,C62,C63,C64,"ATLAS");


      



  //--------------------------------------//
      AnalysisDataPointers total_results;
      //total_results_push_back(atlasData);
      total_results.push_back(cmsData);

      result = total_results;
      
      
    };
    
    void InterpolatedMCInfo(MCLoopInfo& result)
    {
      
      // This makes an MCLoopInfo object for satisfying the LHC
      // likelihood calculation dependency

      result.reset_flags();
      
    }
	

    
  }
    
  }
