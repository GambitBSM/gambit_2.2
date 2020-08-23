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
///  Analyses based on: arxiv:1711.03301 and https://journals.aps.org/prd/abstract/10.1103/PhysRevD.97.092005

//////  *********************************************
#include <chrono>
#include <thread>
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
#include "gambit/Utils/file_lock.hpp"
//#include DMEFT_grids.hpp

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
    // Forward declaration of funtion in LHC_likelihoods
    AnalysisLogLikes calc_loglikes_for_analysis(const AnalysisData&, bool, bool, bool);


    // ---------------------------------- INTERPOLATION FUNCTIONS ------------------------------------------------
      const char* colliderbitdata_path = GAMBIT_DIR "/ColliderBit/data/"; 
      #define PI 3.14159265
      // Initialize all data
      // static const size_t data_INC           = 15;
      // static const size_t data_SIZE          = pow(data_INC,2);
      // static const size_t data_INC_low       = 4;
      // static const size_t data_SIZE_low      = pow(data_INC_low,2);
      // static const size_t data_INC_d7        = 19;
      // static const size_t data_SIZE_d7       = data_INC_d7;


      // static const size_t cms_bin_size       = 22;
      // static const size_t atlas_bin_size     = 10;

      #define data_INC        15
      #define data_SIZE       225
      #define data_INC_low    4
      #define data_SIZE_low   16
      #define data_INC_d7     19
      #define data_SIZE_d7     19
      #define cms_bin_size    22
      #define atlas_bin_size  10


      const char* met_ATLAS_23               = GAMBIT_DIR "/ColliderBit/data/DMEFT/met_hist_ATLAS_C62_C63.txt";
      const char* met_ATLAS_14               = GAMBIT_DIR "/ColliderBit/data/DMEFT/met_hist_ATLAS_C61_C64.txt";
      const char* met_CMS_23                 = GAMBIT_DIR "/ColliderBit/data/DMEFT/met_hist_CMS_C62_C63.txt";
      const char* met_CMS_14                 = GAMBIT_DIR "/ColliderBit/data/DMEFT/met_hist_CMS_C61_C64.txt";
          
      const char* met_ATLAS_23_low           = GAMBIT_DIR "/ColliderBit/data/DMEFT/met_hist_ATLAS_C62_C63_low.txt";
      const char* met_ATLAS_14_low           = GAMBIT_DIR "/ColliderBit/data/DMEFT/met_hist_ATLAS_C61_C64_low.txt";
      const char* met_CMS_23_low             = GAMBIT_DIR "/ColliderBit/data/DMEFT/met_hist_CMS_C62_C63_low.txt";
      const char* met_CMS_14_low             = GAMBIT_DIR "/ColliderBit/data/DMEFT/met_hist_CMS_C61_C64_low.txt";
          
      const char* met_ATLAS_71               = GAMBIT_DIR "/ColliderBit/data/DMEFT/met_hist_ATLAS_C71.txt";
      const char* met_ATLAS_72               = GAMBIT_DIR "/ColliderBit/data/DMEFT/met_hist_ATLAS_C72.txt";
      const char* met_ATLAS_73               = GAMBIT_DIR "/ColliderBit/data/DMEFT/met_hist_ATLAS_C73.txt";
      const char* met_ATLAS_74               = GAMBIT_DIR "/ColliderBit/data/DMEFT/met_hist_ATLAS_C74.txt";
      const char* met_CMS_71                 = GAMBIT_DIR "/ColliderBit/data/DMEFT/met_hist_CMS_C71.txt";
      const char* met_CMS_72                 = GAMBIT_DIR "/ColliderBit/data/DMEFT/met_hist_CMS_C72.txt";
      const char* met_CMS_73                 = GAMBIT_DIR "/ColliderBit/data/DMEFT/met_hist_CMS_C73.txt";
      const char* met_CMS_74                 = GAMBIT_DIR "/ColliderBit/data/DMEFT/met_hist_CMS_C74.txt";
      // ----------------------------------//     
      // Atlas/CMS analysis arrays 
      // ----------------------------------//

      static const double CMS_OBSNUM[cms_bin_size] = {
                              136865, 74340, 42540, 25316, 15653, 10092, 8298, 4906, 2987, 2032, 1514,
                              926, 557, 316, 233, 172, 101, 65, 46, 26, 31, 29};
      static const double CMS_BKGNUM[cms_bin_size] = {
                                          134500, 73400, 42320, 25490, 15430, 10160, 8480, 4865, 2970, 1915, 1506,
                                          844, 526, 325, 223, 169, 107, 88.1, 52.8, 25.0, 25.5, 26.9
                                             };
      static const double CMS_BKGERR[cms_bin_size] = { 
                                            3700, 2000, 810, 490, 310, 170, 140, 95, 49, 33, 32, 18, 14, 12, 9, 8, 6, 5.3, 3.9, 2.5, 2.6, 2.8
                                                };

      static Eigen::MatrixXd m_BKGCOV(22,22);

      static const double ATLAS_OBSNUM[atlas_bin_size] = {111203,67475,35285,27843,8583,2975,1142,512,223,245};
      static const double ATLAS_BKGNUM[atlas_bin_size] = {111100,67100,33820,27640,8360,2825,1094,463,213,226};
      static const double ATLAS_BKGERR[atlas_bin_size] = {2300  ,1400 ,940  ,610  ,190 ,78  ,33  ,19 ,9  ,16 };

      // Define ATLAS and CMS exclusive signal regions: These are arrays of the MIN met in the bin. 
      static const double METMINS_ATLAS[atlas_bin_size]     = {250., 300., 350., 400., 500., 600., 700., 800., 900., 1000.};
      static const double METMINS_CMS[cms_bin_size]       = {250., 280., 310., 340., 370., 400., 430.,470.,510.,550.,590.,640.,690.,740.,790.,840.,900.,960.,1020.,1090.,1160.,1250.};

    double LinearInterpolation(double y2, double y1, double y, double q1,double q2)
    {
      return  ( 1.0 / (y2-y1) ) * ( (y2 - y)*q1 + (y-y1)*q2 );
    }

    double BilinearInterpolation(double q11, double q12, double q21, double q22, 
                                 double x1, double x2, double y1, double y2, 
                                 double x, double y, double yalpha=0,bool debug=false)
    {
      if (q11 < 0 ){
        q11 = LinearInterpolation(y2,yalpha,y1,-1*q11,q12);
      }

      if (q21 < 0 ){
        q21 = LinearInterpolation(y2,yalpha,y1,-1*q21,q22);
      }

      if (q22 < 0 ){
        q22 = LinearInterpolation(yalpha,y1,y2,q21,-1*q22);
      }

      if (q12 < 0 ){
        q12 = LinearInterpolation(yalpha,y1,y2,q11,-1*q12);
      }

      double x2x1, y2y1, x2x, y2y, yy1, xx1;
      x2x1 = x2 - x1;
      y2y1 = y2 - y1;
      x2x = x2 - x;
      y2y = y2 - y;
      yy1 = y - y1;
      xx1 = x - x1;

      if (debug){
        cout << " Oh dear..." << x2x1 << " "<< y2y1<< " "<< x2x << " " << y2y << " "<< yy1<< " "<<xx1<< " x and y "<< x<< " "<< y<<endl;
        cout << " Oh dear..." << y1 << " "<< y2 << endl;

      }
      return float(1.0) / float(x2x1 * y2y1) * ( 
        q11 * x2x * y2y +
        q21 * xx1 * y2y +
        q12 * x2x * yy1 +
        q22 * xx1 * yy1
      );
    }


    // A function to modify the DMEFT LHC signal prediction for ETmiss bins where ETmiss > Lambda
    void signal_modifier_function(AnalysisData& adata, float lambda, float a)
    {
      int met_bin_size;
      const double* METMINS;

      // Choose experiment
      if (adata.analysis_name.find("ATLAS") != string::npos)
      {
        bool is_ATLAS = true;
        METMINS = METMINS_ATLAS;
        met_bin_size = atlas_bin_size;
      }
      else if (adata.analysis_name.find("CMS") != string::npos)
      {
        bool is_CMS = false;
        METMINS = METMINS_CMS;
        met_bin_size = cms_bin_size;
      }
      else
      {
        ColliderBit_error().raise(LOCAL_INFO, "Unknown analysis encountered in signal_modifier_function!");
      }

      // Modify signals
      for (int bin_index = 0; bin_index < met_bin_size; bin_index++ ) 
      {
        double MET_min = METMINS_ATLAS[bin_index];
        double weight = pow(MET_min / lambda, -a);

        if (lambda < MET_min)
        {
          SignalRegionData& srdata = adata[bin_index];
          srdata.n_signal *= weight;
          srdata.n_signal_at_lumi *= weight;
        }
      } 
    }


    // ---------------------------------------------------- //
    //  Calculate Yields // 
    // ---------------------------------------------------- //  
 
    void Acceptance_CS(double * accep, float m,float O1,float O2, float lambda ,const char* pair, const char* experiment)
    {

      static bool first     = true;
      static bool first_low = true;
      Utils::FileLock mylock("Get_data_once");

      if (m>150){

        static double MET_HIST_CMS_14[data_SIZE][cms_bin_size];
        static double MET_HIST_ATLAS_14[data_SIZE][atlas_bin_size];
        static double MET_HIST_CMS_23[data_SIZE][cms_bin_size];
        static double MET_HIST_ATLAS_23[data_SIZE][atlas_bin_size];
        static double THETA_CMS_14[data_SIZE];
        static double THETA_CMS_23[data_SIZE];
        static double THETA_ATLAS_14[data_SIZE];  
        static double THETA_ATLAS_23[data_SIZE];  
        static double MASS_CMS_14[data_SIZE];
        static double MASS_CMS_23[data_SIZE];
        static double MASS_ATLAS_14[data_SIZE];  
        static double MASS_ATLAS_23[data_SIZE];  
        static double CS_CMS_14[data_SIZE];
        static double CS_CMS_23[data_SIZE];
        static double CS_ATLAS_14[data_SIZE];  
        static double CS_ATLAS_23[data_SIZE];  
        static double nJets[data_SIZE];
        // ----------------------------------//
        // Define just mass and angle arrays // 
        static double theta[data_INC];
        static double mass[data_INC]; 

        
    
        if (first)
        {

          mylock.get_lock();

          cout << "Reading in grids. [Only happens on first itteration per MPI process]."<<endl;
          float var1,var2;
          FILE * fp = fopen(GAMBIT_DIR "/ColliderBit/data/DMEFT/X_Y_ATLAS_C62_C63.txt","r");   // The masses and thetas are the same for each! 
          for (int ll = 0; ll < data_INC; ++ll){
            fscanf(fp,"%f %f", &var1, &var2);
            mass[ll] = var1;
            theta[ll]= var2;
            // cout << mass[ll]<<endl;
          }
          fclose(fp); 
          

          ifstream mb(met_ATLAS_23);

          for(int row = 0; row < data_SIZE; row++) {  
            for(int column = 0; column < atlas_bin_size; column++){
              mb >> MET_HIST_ATLAS_23[row][column];
            }
          }
          mb.close();


          ifstream mba14(met_ATLAS_14);
          for(int row = 0; row < data_SIZE; row++) {  
            for(int column = 0; column < atlas_bin_size; column++){
              mba14 >> MET_HIST_ATLAS_14[row][column];
            }
          }
          mba14.close();


          ifstream mbc23(met_CMS_23);
          for(int row = 0; row < data_SIZE; row++) {  
            for(int column = 0; column < cms_bin_size; column++){
              mbc23 >> MET_HIST_CMS_23[row][column];
            }
          }
          mbc23.close();     


          ifstream mbc14(met_CMS_14);
          for(int row = 0; row < data_SIZE; row++) {  
            for(int column = 0; column < cms_bin_size; column++){
              mbc14 >> MET_HIST_CMS_14[row][column];
            }
          }
          mbc14.close();


          float p1,p2,p3,p4;
          FILE * pp = fopen(GAMBIT_DIR "/ColliderBit/data/DMEFT/grid_output_ATLAS_C61_C64.txt","r");              // file containing numbers in 5 columns 
          for (int ll = 0; ll < data_SIZE; ++ll){
            fscanf(pp,"%f %f %f %f", &p1,&p2,&p3,&p4);
            MASS_ATLAS_14[ll] = p1; 
            THETA_ATLAS_14[ll]= p2;
            nJets[ll]         = p3;
            CS_ATLAS_14[ll]   = p4;   
          }
          fclose(pp);


          float a1,a2,a3,a4;
          FILE * ap = fopen(GAMBIT_DIR "/ColliderBit/data/DMEFT/grid_output_ATLAS_C62_C63.txt","r");             // file containing numbers in 5 columns 
          for (int ll = 0; ll < data_SIZE; ++ll){
            fscanf(ap,"%f %f %f %f", &a1,&a2,&a3,&a4);
            MASS_ATLAS_23[ll] = a1; 
            THETA_ATLAS_23[ll]= a2;
            nJets[ll]         = a3;
            CS_ATLAS_23[ll]   = a4;   
          }
          fclose(ap);
          
          float d1,d2,d3,d4;
          FILE *dp=fopen(GAMBIT_DIR "/ColliderBit/data/DMEFT/grid_output_CMS_C61_C64.txt","r");              // file containing numbers in 5 columns 
          for (int ll = 0; ll < data_SIZE; ++ll){
            fscanf(dp,"%f %f %f %f", &d1,&d2,&d3,&d4);
            MASS_CMS_14[ll] = d1; 
            THETA_CMS_14[ll]= d2;
            nJets[ll]       = d3;
            CS_CMS_14[ll]   = d4;   
          }
          fclose(dp);

          float b1,b2,b3,b4;
          FILE * bp = fopen(GAMBIT_DIR "/ColliderBit/data/DMEFT/grid_output_CMS_C62_C63.txt","r");              // file containing numbers in 5 columns 
          for (int ll = 0; ll < data_SIZE; ++ll){
            fscanf(bp,"%f %f %f %f", &b1,&b2,&b3,&b4);
            MASS_CMS_23[ll] = b1; 
            THETA_CMS_23[ll]= b2;
            nJets[ll]       = b3;
            CS_CMS_23[ll]   = b4;   
          }
          fclose(bp);
          
          mylock.release_lock();

          first = false;

        }



        // Define temp. arrays for storing yields. 
        // cout << "Check things "<<mass[0]<<endl;  
        int met_bin_size;
        double ** MET_HIST = new double*[data_SIZE];
        double THETA[data_SIZE];
        double MASS[data_SIZE];
        double CS[data_SIZE];
        
        // cout << "Check things 2 <<mass[0]"<<endl;  

        if (strcmp(experiment,"ATLAS") == 0 && strcmp(pair,"23") == 0)
        {
          met_bin_size = atlas_bin_size;

          // double** MET_HIST = new double*[data_SIZE];
          for(int i = 0; i < data_SIZE; ++i){
            MET_HIST[i] = new double[met_bin_size];
            MASS[i]     = MASS_ATLAS_23[i];
            THETA[i]    = THETA_ATLAS_23[i];
            CS[i]       = CS_ATLAS_23[i];
          }
          // Assign met histogram to current experiment
          for (int kk = 0; kk<data_SIZE;++kk){
            for (int j = 0; j<met_bin_size;++j){
              MET_HIST[kk][j] = MET_HIST_ATLAS_23[kk][j];
            }
          }
        }

        else if (strcmp(experiment,"CMS") == 0 && strcmp(pair,"23") == 0){
          // std::cout << "BITE" << std::endl;

          met_bin_size = cms_bin_size;

          // double** MET_HIST = new double*[data_SIZE];
          for(int i = 0; i < data_SIZE; ++i){
            MET_HIST[i] = new double[met_bin_size];
            MASS[i]     = MASS_CMS_23[i];
            THETA[i]    = THETA_CMS_23[i];
            CS[i]       = CS_CMS_23[i];
          }
          // Assign met histogram to current experiment
          for (int kk = 0; kk<data_SIZE;++kk){
            for (int j = 0; j<met_bin_size;++j){
              MET_HIST[kk][j] = MET_HIST_CMS_23[kk][j];

            }
          }
        }

        else if (strcmp(experiment,"ATLAS") == 0 && strcmp(pair,"14") == 0){
          met_bin_size = atlas_bin_size;

          // double** MET_HIST = new double*[data_SIZE];
          for(int i = 0; i < data_SIZE; ++i){
            MET_HIST[i] = new double[met_bin_size];
            MASS[i]     = MASS_ATLAS_14[i];
            THETA[i]    = THETA_ATLAS_14[i];
            CS[i]       = CS_ATLAS_14[i];
          }
          // Assign met histogram to current experiment
          for (int kk = 0; kk<data_SIZE;++kk){
            for (int j = 0; j<met_bin_size;++j){
              MET_HIST[kk][j] = MET_HIST_ATLAS_14[kk][j];

            }
          }
        }


        else if (strcmp(experiment,"CMS") == 0 && strcmp(pair,"14") == 0){
          met_bin_size = cms_bin_size;

          // double** MET_HIST = new double*[data_SIZE];
          for(int i = 0; i < data_SIZE; ++i){
            MET_HIST[i] = new double[met_bin_size];
            MASS[i]     = MASS_CMS_14[i];
            THETA[i]    = THETA_CMS_14[i];
            CS[i]       = CS_CMS_14[i];
          }
          // Assign met histogram to current experiment
          for (int kk = 0; kk<data_SIZE;++kk){
            for (int j = 0; j<met_bin_size;++j){
              MET_HIST[kk][j] = MET_HIST_CMS_14[kk][j];

            }
          }
        }
        // cout << "Check things 5"<<mass[0]<<endl;  

          
        // Calculate normalisation
        double Norm,th;

        if (O1==0){
          Norm = pow(O2,2);
          th   = 0;
          // cout << " O1 is zero"<< endl;
        }
        else if (O2==0){
          Norm = pow(O1,2);
          th   = float(PI)/float(2);
          // cout << " O2 is zero"<< endl;
        }
        else{
          th    = 0.5*asin(float(2*O1*O2)/float((pow(O1,2)+pow(O2,2))));

          if (O1*O2 < 0){
            th = th + float(PI);
          }
          // cout << "Theta = "<< th<<endl;
          Norm  = 2*O1*O2/(sin(2.0*th));
        }


        // Checks to go ahead with interpolation
        // cout << "Check things 6"<<mass[0]<<endl;  


        if (m<mass[0] || m>mass[data_INC-1]){
          cout<<" Error! Mass param out of range with value "<< m << " Itterator location = " << mass[data_INC-1]<< " Exiting..."<<endl;
          std::exit(EXIT_SUCCESS);
        }
        else if (th<theta[0] || th>theta[data_INC-1]){
          cout<<" Error! Theta param out of range with value " << th << " Exiting..."<<endl;
          std::exit(EXIT_SUCCESS);
        }
        // cout << " Acceptance_CS DEBUG: 4" << endl;

        
        // Get x1,2 y1,2 : Mass and theta coordinates for interpolation
        double x1,x2,y1,y2;
        int xi,yj;
        for(int ii = 0; ii < data_INC-1; ++ii) {
          if (m >= mass[ii] && m <= mass[ii+1]){
            x1 = mass[ii];
            x2 = mass[ii+1]; 
            xi = ii;
            break;
          }
        }
        for(int jj = 0; jj < data_INC-1; ++jj) {
          if (th >= theta[jj] && th <= theta[jj+1]){
            y1 = theta[jj];
            y2 = theta[jj+1];
            yj = jj;
            break;
          }
        }

        // cout << "Check things 7"<<mass[0]<<endl;  

        // Get C's
        double C11=0.0 ,C12=0.0,C21=0.0,C22=0.0,yalpha=0;

        // Define Q's as array: One Q type for each met bin.

        double Q11[met_bin_size];
        double Q12[met_bin_size];
        double Q21[met_bin_size];
        double Q22[met_bin_size];

        // cout << "Check things 8"<<mass[0]<<endl;  

        // NJets and Cross-section
        // // !!!!!!!!!!!!!!!!!!!! HERE AGAIN BUGS !!!!! << endl;

        // std::cout << "met_bin_size: " << met_bin_size << std::endl;

        for (int Emiss = 0; Emiss < met_bin_size; Emiss++ ) 
        {
          Q11[Emiss] = 0.0;
          Q12[Emiss] = 0.0;
          Q21[Emiss] = 0.0;
          Q22[Emiss] = 0.0;
          
          // cout << " Emiss = "<< Emiss<< " Inital Q's: "<< Q11[Emiss]<<" "<< Q12[Emiss]<<" "<< Q12[Emiss] <<" "<< Q22[Emiss]<<endl;
          while (Q11[Emiss]==0.0 || Q12[Emiss]==0.0 || Q21[Emiss]== 0.0 || Q22[Emiss]==0.0 || C11==0.0 || C12==0.0 || C21== 0.0 || C22==0.0)
          { 
            // cout << Q11[Emiss]<<" "<< Q12[Emiss]<<" "<< Q21[Emiss] <<" "<< Q22[Emiss]<<endl;
            // cout << " X1 Y1 X2 Y2  = " << x1<< "  " << y1<< " " <<x2<< " " << y2 << endl;
            for(int kk = 0; kk < data_SIZE; ++kk) 
            {
              // cout << MASS[kk]<<" "<< THETA[kk]<< " Emiss = "<< Emiss <<"|  |"<<MET_HIST[kk][Emiss]<<" " << kk<< " |     |" << x2<<" " << y2<<" "<< Q11[Emiss]<<" "<< Q12[Emiss]<<" "<< Q12[Emiss] <<" "<< Q22[Emiss]<<endl;
              
              if (MASS[kk]==x1 && THETA[kk]==y1)
              {
                // Q11[Emiss] = nJets[kk];
                  if (MET_HIST[kk][Emiss] < 0){
                    Q11[Emiss] = -1*MET_HIST[kk-1][Emiss];
                    if ( std::isnan(Q11[Emiss])){
                      cout << "NAN in dodgey!!! Emiss = "<< Emiss<< ", " << " X1 Y1 X2 Y2  = " << x1<< "  " << y1<< " " <<x2<< " " << y2 << endl;
                    }
                    C11        = -1*CS[kk-1];
                    yalpha     = THETA[kk-1];
                    // cout << "Have made the hack" << endl;
                  }
                  
                  else {
                    Q11[Emiss] = MET_HIST[kk][Emiss];
                    if ( std::isnan(Q11[Emiss])){
                      cout << "NAN!!! Emiss = "<< Emiss<< ", " << " X1 Y1 X2 Y2  = " << x1<< "  " << y1<< " " <<x2<< " " << y2 << endl;
                    }
                    C11 = CS[kk];
                    // cout << "Q11 = " << Q22[Emiss] << " mass, th = "    << MASS[kk]<< "  "<< THETA[kk]<<endl;
                  } 
                  if (Q11[Emiss] == 0){
                    cout << "Q11 not set" << Q11[Emiss]<< " "<< Emiss<< " X1 Y1 X2 Y2  = " << x1<< "  " << y1<< " " <<x2<< " " << y2 << endl;
                }
              }
              else if (MASS[kk]==x1 && THETA[kk]==y2)
              {
                // cout << "Here in loop. K = "<< kk<< " x1 = "<< x1<< " y2 = "<< y2<< " met_hist = " << MET_HIST[kk][Emiss]<<endl;
                
                if (MET_HIST[kk][Emiss] < 0)
                {
                    Q12[Emiss] = -1*MET_HIST[kk+1][Emiss];
                    if ( std::isnan(Q12[Emiss]))
                    {
                      cout << "NAN in dodgey!!! Emiss = "<< Emiss<< ", " << " X1 Y1 X2 Y2  = " << x1<< "  " << y1<< " " <<x2<< " " << y2 << endl;
                    }

                    C12        = -1*CS[kk+1];
                    yalpha     = THETA[kk+1];
                    // cout << "Have made the hack" << endl;
                }
                else {
                    Q12[Emiss] = MET_HIST[kk][Emiss];
                    if ( std::isnan(Q12[Emiss])){
                      cout << "NAN!!! Emiss = "<< Emiss<< ", " << " X1 Y1 X2 Y2  = " << x1<< "  " << y1<< " " <<x2<< " " << y2 << endl;
                    }

                    C12 = CS[kk];
                }
                if (Q12[Emiss] == 0){
                  cout << "Q12 not set" << Q12[Emiss]<< " "<< Emiss<< " X1 Y1 X2 Y2  = " << x1<< "  " << y1<< " " <<x2<< " " << y2 << endl;
                }
              }

              else if (MASS[kk]==x2 && THETA[kk]==y1)
              {

                if (MET_HIST[kk][Emiss] < 0)
                {
                  Q21[Emiss] = -1*MET_HIST[kk-1][Emiss];

                  if ( std::isnan(Q21[Emiss])){
                    cout << "NAN in dodgey!!! Emiss = "<< Emiss<< ", " << " X1 Y1 X2 Y2  = " << x1<< "  " << y1<< " " <<x2<< " " << y2 << endl;
                  }
                  C21        = -1*CS[kk-1];
                  yalpha     = THETA[kk-1];
                  // cout << "Have made the hack" << endl;
                } 

                else
                {
                  Q21[Emiss] = MET_HIST[kk][Emiss];
                    if ( std::isnan(Q21[Emiss])){
                      cout << "NAN!!! Emiss = "<< Emiss<< ", " << " X1 Y X2 Y2  = " << x1<< "  " << y1<< " " <<x2<< " " << y2 << endl;
                    }

                  C21 = CS[kk];
                }
                if (Q21[Emiss] == 0){
                  cout << "Q21 not set" << Q21[Emiss]<< " "<< Emiss<< " X1 Y1 X2 Y2  = " << x1<< "  " << y1<< " " <<x2<< " " << y2 << endl;
                }
              }  

              else if (MASS[kk]==x2 && THETA[kk]==y2)
              {

                  if (MET_HIST[kk][Emiss] < 0)
                  {
                    Q22[Emiss] = -1*MET_HIST[kk+1][Emiss];

                    if ( std::isnan(Q22[Emiss])){
                      cout << "NAN in dodgey!!! Emiss = "<< Emiss<< ", " << " X1 Y1 X2 Y2  = " << x1<< "  " << y1<< " " <<x2<< " " << y2 << endl;
                    }  
                    C22        = -1*CS[kk+1];
                    yalpha     = THETA[kk+1];
                    // cout << "Have made the hack " << Q22[Emiss]<< " "<< C22<< " "<< yalpha <<  endl;
                  } 

                  else 
                  {
                    Q22[Emiss] = MET_HIST[kk][Emiss];
                      if ( std::isnan(Q22[Emiss])){
                        cout << "NAN!!! Emiss = "<< Emiss<< ", " << " X1 Y1 X2 Y2  = " << x1<< "  " << y1<< " " <<x2<< " " << y2 << endl;
                      }
                    C22 = CS[kk];
                  }
                  // cout << " Q22"<< " "<< Q22[Emiss]<<endl;
                  if (Q22[Emiss] == 0){
                    cout << "Q22 not set" << Q22[Emiss]<< " "<< Emiss<< " X1 Y1 X2 Y2  = " << x1<< "  " << y1<< " " <<x2<< " " << y2 << endl;
                  }
                } 
              } // Loop over kk 
            } // While loop

            // cout << "Exited while loop..." << endl;
            // cout << "Check things 9"<<mass[0]<<endl;  

            // cout << " Acceptance_CS DEBUG: 5 - Fixed" << endl;

            // Luminoscity scaling gets applied at the end...
            double A   = BilinearInterpolation(Q11[Emiss], Q12[Emiss], Q21[Emiss], Q22[Emiss], x1, x2, y1, y2, m, th,yalpha);
            double B   = BilinearInterpolation(C11, C12, C21, C22, x1, x2, y1, y2, m, th,yalpha);
            // double res =  36000.0*float(Norm)*A*float(Norm)*B; 
            double res =  36000.0*A*float(Norm)*B; 
            
            // double res =  Norm*BilinearInterpolation(Q11[Emiss], Q12[Emiss], Q21[Emiss], Q22[Emiss], x1, x2, y1, y2, m, th)*Norm*BilinearInterpolation(C11, C12, C21, C22, x1, x2, y1, y2, m, th); 


            if (std::isnan(res))
            {
              cout << " Test within function: Experiment =  "<< experiment << " res =  "<< res << " Pair  = " << pair <<" CS = "<<Norm*BilinearInterpolation(C11, C12, C21, C22, x1, x2, y1, y2, m, th,0,true)<< " Yield = "<< Norm*BilinearInterpolation(Q11[Emiss], Q12[Emiss], Q21[Emiss], Q22[Emiss], x1, x2, y1, y2, m, th,true) <<" Emiss = "<< Emiss << " Q's: "<< Q11[Emiss]<<" " << Q12[Emiss]<<" " << Q21[Emiss]<<" " <<Q22[Emiss]<<" "<< endl;
              cout << " X1 Y1 X2 Y2  = " << x1<< "  " << y1<< " " <<x2<< " " << y2 << endl;
              cout << "th is the problem = "<< th<<endl;
            }
            // cout << "Check things 10"<<mass[0]<<endl;  
            
            //  cout << "Res = "<< res << " Mass, theta = "<< m <<" , "<<th<<" A = "<<A<<" B = "<<B<<endl;
            
            double lambda_scaling = float(pow(1000.0,4))/float(pow(lambda,4));

            accep[Emiss] = res*lambda_scaling;

          } // Loop over Emiss

          // cout << &accep << std::endl;
          // cout << sizeof(accep) << std::endl;
          // cout << "Check things after accep"<<endl;  

          // for(int j=0; j<22; j++) {
          //   cout << &accep[j] << std::endl;
          //    cout << accep[j] << std::endl;
          // }


          // cout << "Check things b4 filling"<<endl;  

          std::fill_n(THETA,data_SIZE,0);
          std::fill_n(MASS,data_SIZE,0);
          std::fill_n(CS,data_SIZE,0);


          // MET_HIST[data_SIZE][met_bin_size] = {};
              //Free each sub-array
          for(int i = 0; i < data_SIZE; ++i) {
              delete[] MET_HIST[i];   
          }
          //Free the array of pointers
          delete[] MET_HIST;


          // MET_HIST_ATLAS[data_SIZE][atlas_bin_size] = {};
          // MET_HIST_CMS[data_SIZE][cms_bin_size] = {};
          // MET_HIST_ATLAS[data_SIZE][atlas_bin_size] = {};
          // cout << "Check things after fillinge "<<endl;  
      }

      else if (m<=150){
        static double MET_HIST_CMS_14_low[data_SIZE_low][cms_bin_size];
        static double MET_HIST_ATLAS_14_low[data_SIZE_low][atlas_bin_size];
        static double MET_HIST_CMS_23_low[data_SIZE_low][cms_bin_size];
        static double MET_HIST_ATLAS_23_low[data_SIZE_low][atlas_bin_size];
        static double THETA_CMS_14_low[data_SIZE_low];
        static double THETA_CMS_23_low[data_SIZE_low];
        static double THETA_ATLAS_14_low[data_SIZE_low];  
        static double THETA_ATLAS_23_low[data_SIZE_low];  
        static double MASS_CMS_14_low[data_SIZE_low];
        static double MASS_CMS_23_low[data_SIZE_low]; 
        static double MASS_ATLAS_14_low[data_SIZE_low];  
        static double MASS_ATLAS_23_low[data_SIZE_low];  
        static double CS_CMS_14_low[data_SIZE_low];
        static double CS_CMS_23_low[data_SIZE];
        static double CS_ATLAS_14_low[data_SIZE_low];  
        static double CS_ATLAS_23_low[data_SIZE_low];  
        static double nJets_low[data_SIZE_low];
        static double theta_low[data_INC_low];
        static double mass_low[data_INC_low]; 

        
    
        if (first_low)
        {
          mylock.get_lock();

          cout << "Reading in grids. [Only happens on first itteration per MPI process]."<<endl;
          float var1,var2;
          FILE * fp = fopen(GAMBIT_DIR "/ColliderBit/data/DMEFT/X_Y_ATLAS_C62_C63_low.txt","r");   // The masses and thetas are the same for each! 
          for (int ll = 0; ll < data_INC_low; ++ll){
            fscanf(fp,"%f %f", &var1, &var2);
            mass_low[ll] = var1;
            theta_low[ll]= var2;
            // cout << mass[ll]<<endl;
          }
          fclose(fp); 
          

          ifstream mb(met_ATLAS_23_low);

          for(int row = 0; row <data_SIZE_low; row++) {  
            for(int column = 0; column < atlas_bin_size; column++){
              mb >> MET_HIST_ATLAS_23_low[row][column];
            }
          }
          mb.close();


          ifstream mba14(met_ATLAS_14_low);
          for(int row = 0; row <data_SIZE_low; row++) {  
            for(int column = 0; column < atlas_bin_size; column++){
              mba14 >> MET_HIST_ATLAS_14_low[row][column];
            }
          }
          mba14.close();


          ifstream mbc23(met_CMS_23_low);
          for(int row = 0; row <data_SIZE_low; row++) {  
            for(int column = 0; column < cms_bin_size; column++){
              mbc23 >> MET_HIST_CMS_23_low[row][column];
            }
          }
          mbc23.close();     


          ifstream mbc14(met_CMS_14_low);
          for(int row = 0; row <data_SIZE_low; row++) {  
            for(int column = 0; column < cms_bin_size; column++){
              mbc14 >> MET_HIST_CMS_14_low[row][column];
            }
          }
          mbc14.close();


          float p1,p2,p3,p4;
          FILE * pp = fopen(GAMBIT_DIR "/ColliderBit/data/DMEFT/grid_output_ATLAS_C61_C64_low.txt","r");              // file containing numbers in 5 columns 
          for (int ll = 0; ll <data_SIZE_low; ++ll){
            fscanf(pp,"%f %f %f %f", &p1,&p2,&p3,&p4);
            MASS_ATLAS_14_low[ll] = p1; 
            THETA_ATLAS_14_low[ll]= p2;
            nJets_low[ll]         = p3;
            CS_ATLAS_14_low[ll]   = p4;   
          }
          fclose(pp);


          float a1,a2,a3,a4;
          FILE * ap = fopen(GAMBIT_DIR "/ColliderBit/data/DMEFT/grid_output_ATLAS_C62_C63_low.txt","r");             // file containing numbers in 5 columns 
          for (int ll = 0; ll <data_SIZE_low; ++ll){
            fscanf(ap,"%f %f %f %f", &a1,&a2,&a3,&a4);
            MASS_ATLAS_23_low[ll] = a1; 
            THETA_ATLAS_23_low[ll]= a2;
            nJets_low[ll]         = a3;
            CS_ATLAS_23_low[ll]   = a4;   
          }
          fclose(ap);
          
          float d1,d2,d3,d4;
          FILE *dp=fopen(GAMBIT_DIR "/ColliderBit/data/DMEFT/grid_output_CMS_C61_C64_low.txt","r");              // file containing numbers in 5 columns 
          for (int ll = 0; ll <data_SIZE_low; ++ll){
            fscanf(dp,"%f %f %f %f", &d1,&d2,&d3,&d4);
            MASS_CMS_14_low[ll] = d1; 
            THETA_CMS_14_low[ll]= d2;
            nJets_low[ll]       = d3;
            CS_CMS_14_low[ll]   = d4;   
          }
          fclose(dp);

          float b1,b2,b3,b4;
          FILE * bp = fopen(GAMBIT_DIR "/ColliderBit/data/DMEFT/grid_output_CMS_C62_C63_low.txt","r");              // file containing numbers in 5 columns 
          for (int ll = 0; ll <data_SIZE_low; ++ll){
            fscanf(bp,"%f %f %f %f", &b1,&b2,&b3,&b4);
            MASS_CMS_23_low[ll] = b1; 
            THETA_CMS_23_low[ll]= b2;
            nJets_low[ll]       = b3;
            CS_CMS_23_low[ll]   = b4;   
          }
          fclose(bp);
          
          mylock.release_lock();

          first_low = false;

        }

        // Define temp. arrays for storing yields. 
        // cout << "Check things "<<mass[0]<<endl;  
        int met_bin_size;
        double ** MET_HIST = new double*[data_SIZE_low];
        double THETA[data_SIZE_low];
        double MASS[data_SIZE_low];
        double CS[data_SIZE_low];
        
        // cout << "Check things 2 <<mass[0]"<<endl;  

        if (strcmp(experiment,"ATLAS") == 0 && strcmp(pair,"23") == 0)
        {
          met_bin_size = atlas_bin_size;

          // double** MET_HIST = new double*[data_SIZE];
          for(int i = 0; i <data_SIZE_low; ++i){
            MET_HIST[i] = new double[met_bin_size];
            MASS[i]     = MASS_ATLAS_23_low[i];
            THETA[i]    = THETA_ATLAS_23_low[i];
            CS[i]       = CS_ATLAS_23_low[i];
          }
          // Assign met histogram to current experiment
          for (int kk = 0; kk<data_SIZE_low;++kk){
            for (int j = 0; j<met_bin_size;++j){
              MET_HIST[kk][j] = MET_HIST_ATLAS_23_low[kk][j];
            }
          }
        }

        else if (strcmp(experiment,"CMS") == 0 && strcmp(pair,"23") == 0){
          // std::cout << "BITE" << std::endl;

          met_bin_size = cms_bin_size;

          // double** MET_HIST = new double*[data_SIZE];
          for(int i = 0; i <data_SIZE_low; ++i){
            MET_HIST[i] = new double[met_bin_size];
            MASS[i]     = MASS_CMS_23_low[i];
            THETA[i]    = THETA_CMS_23_low[i];
            CS[i]       = CS_CMS_23_low[i];
          }
          // Assign met histogram to current experiment
          for (int kk = 0; kk<data_SIZE_low;++kk){
            for (int j = 0; j<met_bin_size;++j){
              MET_HIST[kk][j] = MET_HIST_CMS_23_low[kk][j];

            }
          }
        }

        else if (strcmp(experiment,"ATLAS") == 0 && strcmp(pair,"14") == 0){
          met_bin_size = atlas_bin_size;

          // double** MET_HIST = new double*[data_SIZE];
          for(int i = 0; i <data_SIZE_low; ++i){
            MET_HIST[i] = new double[met_bin_size];
            MASS[i]     = MASS_ATLAS_14_low[i];
            THETA[i]    = THETA_ATLAS_14_low[i];
            CS[i]       = CS_ATLAS_14_low[i];
          }
          // Assign met histogram to current experiment
          for (int kk = 0; kk<data_SIZE_low;++kk){
            for (int j = 0; j<met_bin_size;++j){
              MET_HIST[kk][j] = MET_HIST_ATLAS_14_low[kk][j];

            }
          }
        }


        else if (strcmp(experiment,"CMS") == 0 && strcmp(pair,"14") == 0){
          met_bin_size = cms_bin_size;

          // double** MET_HIST = new double*[data_SIZE];
          for(int i = 0; i <data_SIZE_low; ++i){
            MET_HIST[i] = new double[met_bin_size];
            MASS[i]     = MASS_CMS_14_low[i];
            THETA[i]    = THETA_CMS_14_low[i];
            CS[i]       = CS_CMS_14_low[i];
          }
          // Assign met histogram to current experiment
          for (int kk = 0; kk<data_SIZE_low;++kk){
            for (int j = 0; j<met_bin_size;++j){
              MET_HIST[kk][j] = MET_HIST_CMS_14_low[kk][j];

            }
          }
        }
        // cout << "Check things 5"<<mass[0]<<endl;  

          
        // Calculate normalisation
        double Norm,th;

        if (O1==0){
          Norm = pow(O2,2);
          th   = 0;
          // cout << " O1 is zero"<< endl;
        }
        else if (O2==0){
          Norm = pow(O1,2);
          th   = float(PI)/2;
          // cout << " O2 is zero"<< endl;
        }
        else{
          th    = 0.5*asin(float(2*O1*O2)/float((pow(O1,2)+pow(O2,2))));

          if (O1*O2 < 0){
            th = th + float(PI);
          }
          // cout << "Theta = "<< th<<endl;
          Norm  = 2*O1*O2/(sin(2.0*th));
        }


        // Checks to go ahead with interpolation
        // cout << "Check things 6"<<mass[0]<<endl;  


        if (m<mass_low[0] || m>mass_low[data_INC_low-1]){
          cout<<" Error! Mass param out of range with value "<< m << " Itterator location = " << mass_low[data_INC_low-1]<< " Exiting..."<<endl;
          std::exit(EXIT_SUCCESS);
        }
        else if (th<theta_low[0] || th>theta_low[data_INC_low-1]){
          cout<<" Error! Theta param out of range with value " << th << " Exiting..."<<endl;
          std::exit(EXIT_SUCCESS);
        }
        // cout << " Acceptance_CS DEBUG: 4" << endl;

        
        // Get x1,2 y1,2 : Mass and theta coordinates for interpolation
        double x1,x2,y1,y2;
        int xi,yj;
        for(int ii = 0; ii <data_INC_low-1; ++ii) {
          if (m >= mass_low[ii] && m <= mass_low[ii+1]){
            x1 = mass_low[ii];
            x2 = mass_low[ii+1]; 
            xi = ii;
            break;
          }
        }
        for(int jj = 0; jj <data_INC_low-1; ++jj) {
          if (th >= theta_low[jj] && th <= theta_low[jj+1]){
            y1 = theta_low[jj];
            y2 = theta_low[jj+1];
            yj = jj;
            break;
          }
        }

        // cout << "Check things 7"<<mass[0]<<endl;  

        // Get C's
        double C11=0.0 ,C12=0.0,C21=0.0,C22=0.0,yalpha=0;

        // Define Q's as array: One Q type for each met bin.

        double Q11[met_bin_size];
        double Q12[met_bin_size];
        double Q21[met_bin_size];
        double Q22[met_bin_size];

        // cout << "Check things 8"<<mass[0]<<endl;  

        // NJets and Cross-section
        // // !!!!!!!!!!!!!!!!!!!! HERE AGAIN BUGS !!!!! << endl;

        // std::cout << "met_bin_size: " << met_bin_size << std::endl;

        for (int Emiss = 0; Emiss < met_bin_size; Emiss++ ) 
        {
          Q11[Emiss] = 0.0;
          Q12[Emiss] = 0.0;
          Q21[Emiss] = 0.0;
          Q22[Emiss] = 0.0;
          
          // cout << " Emiss = "<< Emiss<< " Inital Q's: "<< Q11[Emiss]<<" "<< Q12[Emiss]<<" "<< Q12[Emiss] <<" "<< Q22[Emiss]<<endl;
          while (Q11[Emiss]==0.0 || Q12[Emiss]==0.0 || Q21[Emiss]== 0.0 || Q22[Emiss]==0.0 || C11==0.0 || C12==0.0 || C21== 0.0 || C22==0.0)
          { 
            // cout << Q11[Emiss]<<" "<< Q12[Emiss]<<" "<< Q21[Emiss] <<" "<< Q22[Emiss]<<endl;
            // cout << " X1 Y1 X2 Y2  = " << x1<< "  " << y1<< " " <<x2<< " " << y2 << endl;
            for(int kk = 0; kk <data_SIZE_low; ++kk) 
            {
              // cout << MASS[kk]<<" "<< THETA[kk]<< " Emiss = "<< Emiss <<"|  |"<<MET_HIST[kk][Emiss]<<" " << kk<< " |     |" << x2<<" " << y2<<" "<< Q11[Emiss]<<" "<< Q12[Emiss]<<" "<< Q12[Emiss] <<" "<< Q22[Emiss]<<endl;
              
              if (MASS[kk]==x1 && THETA[kk]==y1)
              {
                // Q11[Emiss] = nJets[kk];
                  if (MET_HIST[kk][Emiss] < 0){
                    Q11[Emiss] = -1*MET_HIST[kk-1][Emiss];
                    if ( std::isnan(Q11[Emiss])){
                      cout << "NAN in dodgey!!! Emiss = "<< Emiss<< ", " << " X1 Y1 X2 Y2  = " << x1<< "  " << y1<< " " <<x2<< " " << y2 << endl;
                    }
                    C11        = -1*CS[kk-1];
                    yalpha     = THETA[kk-1];
                    // cout << "Have made the hack" << endl;
                  }
                  
                  else {
                    Q11[Emiss] = MET_HIST[kk][Emiss];
                    if ( std::isnan(Q11[Emiss])){
                      cout << "NAN!!! Emiss = "<< Emiss<< ", " << " X1 Y1 X2 Y2  = " << x1<< "  " << y1<< " " <<x2<< " " << y2 << endl;
                    }
                    C11 = CS[kk];
                    // cout << "Q11 = " << Q22[Emiss] << " mass, th = "    << MASS[kk]<< "  "<< THETA[kk]<<endl;
                  } 
                  if (Q11[Emiss] == 0){
                    cout << "Q11 not set" << Q11[Emiss]<< " "<< Emiss<< " X1 Y1 X2 Y2  = " << x1<< "  " << y1<< " " <<x2<< " " << y2 << endl;
                }
              }
              else if (MASS[kk]==x1 && THETA[kk]==y2)
              {
                // cout << "Here in loop. K = "<< kk<< " x1 = "<< x1<< " y2 = "<< y2<< " met_hist = " << MET_HIST[kk][Emiss]<<endl;
                
                if (MET_HIST[kk][Emiss] < 0)
                {
                    Q12[Emiss] = -1*MET_HIST[kk+1][Emiss];
                    if ( std::isnan(Q12[Emiss]))
                    {
                      cout << "NAN in dodgey!!! Emiss = "<< Emiss<< ", " << " X1 Y1 X2 Y2  = " << x1<< "  " << y1<< " " <<x2<< " " << y2 << endl;
                    }

                    C12        = -1*CS[kk+1];
                    yalpha     = THETA[kk+1];
                    // cout << "Have made the hack" << endl;
                }
                else {
                    Q12[Emiss] = MET_HIST[kk][Emiss];
                    if ( std::isnan(Q12[Emiss])){
                      cout << "NAN!!! Emiss = "<< Emiss<< ", " << " X1 Y1 X2 Y2  = " << x1<< "  " << y1<< " " <<x2<< " " << y2 << endl;
                    }

                    C12 = CS[kk];
                }
                if (Q12[Emiss] == 0){
                  cout << "Q12 not set" << Q12[Emiss]<< " "<< Emiss<< " X1 Y1 X2 Y2  = " << x1<< "  " << y1<< " " <<x2<< " " << y2 << endl;
                }
              }

              else if (MASS[kk]==x2 && THETA[kk]==y1)
              {

                if (MET_HIST[kk][Emiss] < 0)
                {
                  Q21[Emiss] = -1*MET_HIST[kk-1][Emiss];

                  if ( std::isnan(Q21[Emiss])){
                    cout << "NAN in dodgey!!! Emiss = "<< Emiss<< ", " << " X1 Y1 X2 Y2  = " << x1<< "  " << y1<< " " <<x2<< " " << y2 << endl;
                  }
                  C21        = -1*CS[kk-1];
                  yalpha     = THETA[kk-1];
                  // cout << "Have made the hack" << endl;
                } 

                else
                {
                  Q21[Emiss] = MET_HIST[kk][Emiss];
                    if ( std::isnan(Q21[Emiss])){
                      cout << "NAN!!! Emiss = "<< Emiss<< ", " << " X1 Y X2 Y2  = " << x1<< "  " << y1<< " " <<x2<< " " << y2 << endl;
                    }

                  C21 = CS[kk];
                }
                if (Q21[Emiss] == 0){
                  cout << "Q21 not set" << Q21[Emiss]<< " "<< Emiss<< " X1 Y1 X2 Y2  = " << x1<< "  " << y1<< " " <<x2<< " " << y2 << endl;
                }
              }  

              else if (MASS[kk]==x2 && THETA[kk]==y2)
              {

                  if (MET_HIST[kk][Emiss] < 0)
                  {
                    Q22[Emiss] = -1*MET_HIST[kk+1][Emiss];

                    if ( std::isnan(Q22[Emiss])){
                      cout << "NAN in dodgey!!! Emiss = "<< Emiss<< ", " << " X1 Y1 X2 Y2  = " << x1<< "  " << y1<< " " <<x2<< " " << y2 << endl;
                    }  
                    C22        = -1*CS[kk+1];
                    yalpha     = THETA[kk+1];
                    // cout << "Have made the hack " << Q22[Emiss]<< " "<< C22<< " "<< yalpha <<  endl;
                  } 

                  else 
                  {
                    Q22[Emiss] = MET_HIST[kk][Emiss];
                      if ( std::isnan(Q22[Emiss])){
                        cout << "NAN!!! Emiss = "<< Emiss<< ", " << " X1 Y1 X2 Y2  = " << x1<< "  " << y1<< " " <<x2<< " " << y2 << endl;
                      }
                    C22 = CS[kk];
                  }
                  // cout << " Q22"<< " "<< Q22[Emiss]<<endl;
                  if (Q22[Emiss] == 0){
                    cout << "Q22 not set" << Q22[Emiss]<< " "<< Emiss<< " X1 Y1 X2 Y2  = " << x1<< "  " << y1<< " " <<x2<< " " << y2 << endl;
                  }
                } 
              } // Loop over kk 
            } // While loop

            // cout << "Exited while loop..." << endl;
            // cout << "Check things 9"<<mass[0]<<endl;  

            // cout << " Acceptance_CS DEBUG: 5 - Fixed" << endl;

            // Luminoscity scaling gets applied at the end...
            double A   = BilinearInterpolation(Q11[Emiss], Q12[Emiss], Q21[Emiss], Q22[Emiss], x1, x2, y1, y2, m, th,yalpha);
            double B   = BilinearInterpolation(C11, C12, C21, C22, x1, x2, y1, y2, m, th,yalpha);
            // double res =  36000.0*Norm*A*Norm*B;
            double res =  36000.0*A*float(Norm)*B; 

            // double res =  Norm*BilinearInterpolation(Q11[Emiss], Q12[Emiss], Q21[Emiss], Q22[Emiss], x1, x2, y1, y2, m, th)*Norm*BilinearInterpolation(C11, C12, C21, C22, x1, x2, y1, y2, m, th); 
          
            if (std::isnan(res))
            {
              cout << " Test within function: Experiment =  "<< experiment << " res =  "<< res << " Pair  = " << pair <<" CS = "<<Norm*BilinearInterpolation(C11, C12, C21, C22, x1, x2, y1, y2, m, th,0,true)<< " Yield = "<< Norm*BilinearInterpolation(Q11[Emiss], Q12[Emiss], Q21[Emiss], Q22[Emiss], x1, x2, y1, y2, m, th,true) <<" Emiss = "<< Emiss << " Q's: "<< Q11[Emiss]<<" " << Q12[Emiss]<<" " << Q21[Emiss]<<" " <<Q22[Emiss]<<" "<< endl;
              cout << " X1 Y1 X2 Y2  = " << x1<< "  " << y1<< " " <<x2<< " " << y2 << endl;
              cout << "th is the problem = "<< th<<endl;
            }
            // cout << "Check things 10"<<mass[0]<<endl;  
            
            //  cout << "Res = "<< res << " Mass, theta = "<< m <<" , "<<th<<" A = "<<A<<" B = "<<B<<endl;
            
            double lambda_scaling = float(pow(1000.0,4))/float(pow(lambda,4));

            accep[Emiss] = res*lambda_scaling;
 

          } // Loop over Emiss

          // cout << &accep << std::endl;
          // cout << sizeof(accep) << std::endl;
          // cout << "Check things after accep"<<endl;  

          // for(int j=0; j<22; j++) {
          //   cout << &accep[j] << std::endl;
          //    cout << accep[j] << std::endl;
          // }


          // cout << "Check things b4 filling"<<endl;  

          std::fill_n(THETA,data_SIZE,0);
          std::fill_n(MASS,data_SIZE,0);
          std::fill_n(CS,data_SIZE,0);


          // MET_HIST[data_SIZE][met_bin_size] = {};
              //Free each sub-array
          for(int i = 0; i <data_SIZE_low; ++i) {
              delete[] MET_HIST[i];   
          }
          //Free the array of pointers
          delete[] MET_HIST;

      }
    
    
    
    } // Function Acceptance_CS <---- End of function


    void Acceptance_CS_dim7(double * accep, float m, float Opp, float lambda ,const char* pair, const char* experiment)
    {
      static bool first_7     = true;
      Utils::FileLock mylock7("Get_data_once7");

  
      static double MET_HIST_CMS_71[data_SIZE_d7][cms_bin_size];
      static double MET_HIST_ATLAS_71[data_SIZE_d7][atlas_bin_size];
      static double MET_HIST_CMS_72[data_SIZE_d7][cms_bin_size];
      static double MET_HIST_ATLAS_72[data_SIZE_d7][atlas_bin_size];
      static double MET_HIST_CMS_73[data_SIZE_d7][cms_bin_size];
      static double MET_HIST_ATLAS_73[data_SIZE_d7][atlas_bin_size];
      static double MET_HIST_CMS_74[data_SIZE_d7][cms_bin_size];
      static double MET_HIST_ATLAS_74[data_SIZE_d7][atlas_bin_size];
      static double CS_CMS_71[data_SIZE_d7];
      static double CS_ATLAS_71[data_SIZE_d7];  
      static double CS_CMS_72[data_SIZE_d7];
      static double CS_ATLAS_72[data_SIZE_d7];  
      static double CS_CMS_73[data_SIZE_d7];
      static double CS_ATLAS_73[data_SIZE_d7];  
      static double CS_CMS_74[data_SIZE_d7];
      static double CS_ATLAS_74[data_SIZE_d7];  
      static double mass[data_INC_d7]; 

  
      if (first_7)
      {
        
        mylock7.get_lock();
        
        cout << "Reading in grids. [Only happens on first itteration per MPI process]."<<endl;
        float var1;
        FILE * fp = fopen(GAMBIT_DIR "/ColliderBit/data/DMEFT/X_Y_ATLAS_C71.txt","r");   // The masses and thetas are the same for each! 
        for (int ll = 0; ll < data_INC_d7; ++ll){
          fscanf(fp,"%f", &var1);
          mass[ll] = var1;
        }
        fclose(fp); 
        

        ifstream mb(met_ATLAS_71);
        for(int row = 0; row < data_SIZE_d7; row++) {  
          for(int column = 0; column < atlas_bin_size; column++){
            mb >> MET_HIST_ATLAS_71[row][column];
          }
        }
        mb.close();


        ifstream mba72(met_ATLAS_72);
        for(int row = 0; row < data_SIZE_d7; row++) {  
          for(int column = 0; column < atlas_bin_size; column++){
            mba72 >> MET_HIST_ATLAS_72[row][column];
          }
        }
        mba72.close();

        ifstream mb73(met_ATLAS_73);
        for(int row = 0; row < data_SIZE_d7; row++) {  
          for(int column = 0; column < atlas_bin_size; column++){
            mb73 >> MET_HIST_ATLAS_73[row][column];
          }
        }
        mb73.close();

        ifstream mba74(met_ATLAS_74);
        for(int row = 0; row < data_SIZE_d7; row++) {  
          for(int column = 0; column < atlas_bin_size; column++){
            mba74 >> MET_HIST_ATLAS_74[row][column];
          }
        }
        mba74.close();


        // -------------------------------- //

        ifstream mbc71(met_CMS_71);
        for(int row = 0; row < data_SIZE_d7; row++) {  
          for(int column = 0; column < cms_bin_size; column++){
            mbc71 >> MET_HIST_CMS_71[row][column];
          }
        }
        mbc71.close();


        ifstream mbac72(met_CMS_72);
        for(int row = 0; row < data_SIZE_d7; row++) {  
          for(int column = 0; column < cms_bin_size; column++){
            mbac72 >> MET_HIST_CMS_72[row][column];
          }
        }
        mbac72.close();

        ifstream mbc73(met_CMS_73);
        for(int row = 0; row < data_SIZE_d7; row++) {  
          for(int column = 0; column < cms_bin_size; column++){
            mbc73 >> MET_HIST_CMS_73[row][column];
          }
        }
        mbc73.close();

        ifstream mbac74(met_CMS_74);
        for(int row = 0; row < data_SIZE_d7; row++) {  
          for(int column = 0; column < cms_bin_size; column++){
            mbac74 >> MET_HIST_CMS_74[row][column];
          }
        }
        mbac74.close();

        // Segfault here ^^^^^^^^^^^^^^^^^^^^

        // -------------------------------------------------//
        float p1,p2;
        FILE * pp = fopen(GAMBIT_DIR "/ColliderBit/data/DMEFT/grid_output_ATLAS_C71.txt","r");              // file containing numbers in 5 columns 
        for (int ll = 0; ll < data_SIZE_d7; ++ll){
          fscanf(pp,"%f %f", &p1,&p2);
          CS_ATLAS_71[ll]   = p2;   
        }
        fclose(pp);


        float a1,a4;
        FILE * ap = fopen(GAMBIT_DIR "/ColliderBit/data/DMEFT/grid_output_ATLAS_C72.txt","r");             // file containing numbers in 5 columns 
        for (int ll = 0; ll < data_SIZE_d7; ++ll){
          fscanf(ap,"%f %f", &a1,&a4);
          CS_ATLAS_72[ll]   = a4;   
        }
        fclose(ap);


        float p4;
        pp = fopen(GAMBIT_DIR "/ColliderBit/data/DMEFT/grid_output_ATLAS_C73.txt","r");              // file containing numbers in 5 columns 
        for (int ll = 0; ll < data_SIZE_d7; ++ll){
          fscanf(pp,"%f %f", &p1,&p4);
          CS_ATLAS_73[ll]   = p4;   
        }
        fclose(pp);


        ap = fopen(GAMBIT_DIR "/ColliderBit/data/DMEFT/grid_output_ATLAS_C74.txt","r");             // file containing numbers in 5 columns 
        for (int ll = 0; ll < data_SIZE_d7; ++ll){
          fscanf(ap,"%f %f", &a1,&a4);
          CS_ATLAS_74[ll]   = a4;   
        }
        fclose(ap);
        
       
        pp = fopen(GAMBIT_DIR "/ColliderBit/data/DMEFT/grid_output_CMS_C71.txt","r");              // file containing numbers in 5 columns 
        for (int ll = 0; ll < data_SIZE_d7; ++ll){
          fscanf(pp,"%f %f", &p1,&p2);
          CS_CMS_71[ll]   = p2;   
        }
        fclose(pp);


        ap = fopen(GAMBIT_DIR "/ColliderBit/data/DMEFT/grid_output_CMS_C72.txt","r");             // file containing numbers in 5 columns 
        for (int ll = 0; ll < data_SIZE_d7; ++ll){
          fscanf(ap,"%f %f", &a1,&a4);
          CS_CMS_72[ll]   = a4;   
        }
        fclose(ap);

        pp = fopen(GAMBIT_DIR "/ColliderBit/data/DMEFT/grid_output_CMS_C73.txt","r");              // file containing numbers in 5 columns 
        for (int ll = 0; ll < data_SIZE_d7; ++ll){
          fscanf(pp,"%f %f", &p1,&p4);
          CS_CMS_73[ll]   = p4;   
        }
        fclose(pp);


        ap = fopen(GAMBIT_DIR "/ColliderBit/data/DMEFT/grid_output_CMS_C74.txt","r");             // file containing numbers in 5 columns 
        for (int ll = 0; ll < data_SIZE_d7; ++ll){
          fscanf(ap,"%f %f", &a1,&a4);
          CS_CMS_74[ll]   = a4;   
        }
        fclose(ap);

        mylock7.release_lock();

        first_7 = false;

      }



      int met_bin_size;
      double ** MET_HIST = new double*[data_SIZE_d7];
      double CS[data_SIZE_d7];
      



      if (strcmp(experiment,"ATLAS") == 0 && strcmp(pair,"71") == 0)
      {
        met_bin_size = atlas_bin_size;

        // double** MET_HIST = new double*[data_SIZE_d7];
        for(int i = 0; i < data_SIZE_d7; ++i){
          MET_HIST[i] = new double[met_bin_size];
          CS[i]       = CS_ATLAS_71[i];
        }
        // Assign met histogram to current experiment
        for (int kk = 0; kk<data_SIZE_d7;++kk){
          for (int j = 0; j<met_bin_size;++j){
            MET_HIST[kk][j] = MET_HIST_ATLAS_71[kk][j];
          }
        }
      }

      else if (strcmp(experiment,"CMS") == 0 && strcmp(pair,"71") == 0){
        // std::cout << "BITE" << std::endl;

        met_bin_size = cms_bin_size;

        // double** MET_HIST = new double*[data_SIZE_d7];
        for(int i = 0; i < data_SIZE_d7; ++i){
          MET_HIST[i] = new double[met_bin_size];
          CS[i]       = CS_CMS_71[i];
        }
        // Assign met histogram to current experiment
        for (int kk = 0; kk<data_SIZE_d7;++kk){
          for (int j = 0; j<met_bin_size;++j){
            MET_HIST[kk][j] = MET_HIST_CMS_71[kk][j];

          }
        }
      }

      else if (strcmp(experiment,"ATLAS") == 0 && strcmp(pair,"72") == 0){
        met_bin_size = atlas_bin_size;

        // double** MET_HIST = new double*[data_SIZE_d7];
        for(int i = 0; i < data_SIZE_d7; ++i){
          MET_HIST[i] = new double[met_bin_size];
          CS[i]       = CS_ATLAS_72[i];
        }
        // Assign met histogram to current experiment
        for (int kk = 0; kk<data_SIZE_d7;++kk){
          for (int j = 0; j<met_bin_size;++j){
            MET_HIST[kk][j] = MET_HIST_ATLAS_72[kk][j];

          }
        }
      }


      else if (strcmp(experiment,"CMS") == 0 && strcmp(pair,"72") == 0){
        met_bin_size = cms_bin_size;

        // double** MET_HIST = new double*[data_SIZE_d7];
        for(int i = 0; i < data_SIZE_d7; ++i){
          MET_HIST[i] = new double[met_bin_size];
          CS[i]       = CS_CMS_72[i];
        }
        // Assign met histogram to current experiment
        for (int kk = 0; kk<data_SIZE_d7;++kk){
          for (int j = 0; j<met_bin_size;++j){
            MET_HIST[kk][j] = MET_HIST_CMS_72[kk][j];

          }
        }
      }

      else if (strcmp(experiment,"ATLAS") == 0 && strcmp(pair,"73") == 0){
        met_bin_size = atlas_bin_size;

        // double** MET_HIST = new double*[data_SIZE_d7];
        for(int i = 0; i < data_SIZE_d7; ++i){
          MET_HIST[i] = new double[met_bin_size];
          CS[i]       = CS_ATLAS_73[i];
        }
        // Assign met histogram to current experiment
        for (int kk = 0; kk<data_SIZE_d7;++kk){
          for (int j = 0; j<met_bin_size;++j){
            MET_HIST[kk][j] = MET_HIST_ATLAS_73[kk][j];

          }
        }
      }


      else if (strcmp(experiment,"CMS") == 0 && strcmp(pair,"73") == 0){
        met_bin_size = cms_bin_size;

        // double** MET_HIST = new double*[data_SIZE_d7];
        for(int i = 0; i < data_SIZE_d7; ++i){
          MET_HIST[i] = new double[met_bin_size];
          CS[i]       = CS_CMS_73[i];
        }
        // Assign met histogram to current experiment
        for (int kk = 0; kk<data_SIZE_d7;++kk){
          for (int j = 0; j<met_bin_size;++j){
            MET_HIST[kk][j] = MET_HIST_CMS_73[kk][j];

          }
        }
      }


      else if (strcmp(experiment,"ATLAS") == 0 && strcmp(pair,"74") == 0){
        met_bin_size = atlas_bin_size;

        // double** MET_HIST = new double*[data_SIZE_d7];
        for(int i = 0; i < data_SIZE_d7; ++i){
          MET_HIST[i] = new double[met_bin_size];
          CS[i]       = CS_ATLAS_74[i];
        }
        // Assign met histogram to current experiment
        for (int kk = 0; kk<data_SIZE_d7;++kk){
          for (int j = 0; j<met_bin_size;++j){
            MET_HIST[kk][j] = MET_HIST_ATLAS_74[kk][j];

          }
        }
      }


      else if (strcmp(experiment,"CMS") == 0 && strcmp(pair,"74") == 0){
        met_bin_size = cms_bin_size;

        // double** MET_HIST = new double*[data_SIZE_d7];
        for(int i = 0; i < data_SIZE_d7; ++i){
          MET_HIST[i] = new double[met_bin_size];
          CS[i]       = CS_CMS_74[i];
        }
        // Assign met histogram to current experiment
        for (int kk = 0; kk<data_SIZE_d7;++kk){
          for (int j = 0; j<met_bin_size;++j){
            MET_HIST[kk][j] = MET_HIST_CMS_74[kk][j];

          }
        }
      }

      if (m<mass[0] || m>mass[data_INC_d7-1]){
        cout<<" Error! Mass param out of range with value "<< m << " Itterator location = " << mass[data_INC_d7-1]<< " Exiting..."<<endl;
        std::exit(EXIT_SUCCESS);
      }

      double x1,x2;
      int xi,yj;
      for(int ii = 0; ii < data_INC_d7-1; ++ii) {
        if (m >= mass[ii] && m <= mass[ii+1]){
          x1 = mass[ii];
          x2 = mass[ii+1]; 
          xi = ii;
          break;
        }
      }

      // Define C's
      double C1,C2;
      // Define Q's as array.
      double Q1[met_bin_size];
      double Q2[met_bin_size];

      // Debugging: 
      //int itt = 0;

      for (int Emiss = 0; Emiss < met_bin_size; Emiss++ ) 
      {
        Q1[Emiss] = 0.0;
        Q2[Emiss] = 0.0;
        C1        =0.0;
        C2        =0.0;
        while (Q1[Emiss]==0.0 || Q2[Emiss]==0.0 || C1==0.0 || C2==0.0)

        { 
          // cout<< "Before " << Q1[Emiss]<<" "<< Q2[Emiss]<< " "<<C1<<" "<<C2<<endl;
          // cout << " X1 X2 Emiss metbinsize = " << x1<< " "<<x2<<  "  "<< Emiss<<"  " << met_bin_size<< endl;
          for(int kk = 0; kk < data_SIZE_d7; ++kk) 
          {
              if (mass[kk]==x1)
              {
                // cout << "in x1 ... "<<CS[kk] << " "<<  MET_HIST[kk][Emiss] << endl;
                Q1[Emiss] = MET_HIST[kk][Emiss];
                C1 = CS[kk];
              }
           //----------------- 
              else if (mass[kk]==x2)
              {
                // cout << "in x2 ... "<<CS[kk] << " "<<  MET_HIST[kk][Emiss] << endl;
                  Q2[Emiss] = MET_HIST[kk][Emiss];
                  C2 = CS[kk];
              }  
          } // Loop over kk


          // cout<< "After q's " << Q1[Emiss]<<" "<< Q2[Emiss]<<endl;
          // cout<< "After c's " << C1<<" "<< C2<<endl;
          // cout<< "After Emiss " << Emiss<<endl;

          // itt = itt +1; // Useful for debugging this loop.
          // if (itt==met_bin_size+1){
          //   exit(0);
          // }


        } // While loop






          // Luminoscity scaling gets applied at the end...

          double A   = LinearInterpolation(x2,x1,m,Q1[Emiss],Q2[Emiss]);
          double B   = LinearInterpolation(x2,x1,m,C1,C2);
          double Norm= pow(Opp,2);
          // double res =  36000.0*Norm*A*Norm*B; 
          double res =  36000.0*A*Norm*B; 
          double lambda_scaling = float(pow(1000.0,6))/float(pow(lambda,6));
  
          accep[Emiss] = res*lambda_scaling;

        } // Loop over Emiss


        std::fill_n(CS,data_SIZE_d7,0);

        for(int i = 0; i < data_SIZE_d7; ++i) {
            delete[] MET_HIST[i];   
        }
        //Free the array of pointers
        delete[] MET_HIST;


      }

  
   // Function Acceptance_CS <---- End of function




    void L_Acc_Eff_CS(double * YIELDS, float m,float C61,float C62,float C63, float C64,float C71,float C72,float C73, float C74,float lambda , const char* exper_)
    {
        // cout << "Check things 11 "<<endl;  

      char const *tt = "23";
      char const *of = "14";
      char const *so = "71";
      char const *st = "72";
      char const *sth= "73";
      char const *sf = "74";            

      int met_bin_size;

      if (strcmp(exper_,"ATLAS") == 0){
        met_bin_size = atlas_bin_size;
      }
      else if (strcmp(exper_,"CMS") == 0){
        met_bin_size = cms_bin_size;
      }

      // Dim-6 yields
      double A23[met_bin_size];
      double A14[met_bin_size];

      Acceptance_CS(A23, m, C62, C63, lambda, tt, exper_);

      Acceptance_CS(A14, m, C61, C64, lambda, of, exper_);

      // Dim-7 yields
      double A71[met_bin_size];
      double A72[met_bin_size];
      double A73[met_bin_size];
      double A74[met_bin_size];

      Acceptance_CS_dim7(A71, m, C71, lambda, so, exper_);

      Acceptance_CS_dim7(A72, m, C72, lambda, st, exper_);

      Acceptance_CS_dim7(A73, m, C73, lambda, sth, exper_);

      Acceptance_CS_dim7(A74, m, C74, lambda, sf, exper_);

      // Add yields linearly
      for (int ii = 0; ii < met_bin_size; ++ii){
        YIELDS[ii] = A23[ii] + A14[ii] + A71[ii] + A72[ii] + A73[ii] + A74[ii];
      }

    }

    void DMEFT_results(AnalysisDataPointers& result)
    { 

      // auto start_wall_clock = std::chrono::steady_clock::now();


      using namespace Pipes::DMEFT_results;

      // Clear previous vectors, etc.
      result.clear();

      // Create the thread_local AnalysisData instances we need, 
      // and make sure they are properly cleared for each new point
      thread_local AnalysisData cmsData("CMS_13TeV_MONOJET_36invfb_interpolated");
      cmsData.clear();

      thread_local AnalysisData cmsData_nocovar("CMS_13TeV_MONOJET_36invfb_interpolated_nocovar");
      cmsData_nocovar.clear();

      thread_local AnalysisData atlasData("ATLAS_13TeV_MONOJET_36invfb_interpolated");
      atlasData.clear();


      // cout << "void is run"<< endl;

      const Spectrum& spec = *Dep::DMEFT_spectrum;

      // TODO change floats -> doubles
      float C61 = spec.get(Par::dimensionless, "C61");
      float C62 = spec.get(Par::dimensionless, "C62");
      float C63 = spec.get(Par::dimensionless, "C63");
      float C64 = spec.get(Par::dimensionless, "C64");
      float C71 = spec.get(Par::dimensionless, "C71");
      float C72 = spec.get(Par::dimensionless, "C72");
      float C73 = spec.get(Par::dimensionless, "C73");
      float C74 = spec.get(Par::dimensionless, "C74");
      float lambda = spec.get(Par::mass1, "Lambda");      
      float mchi = spec.get(Par::Pole_Mass, "chi");

      // Do not get segfault when I do a get data here.....
      // float C61   = *Pipes::DMEFT_results::Param["C61"]; 
      // float C62   = *Pipes::DMEFT_results::Param["C62"];
      // float C63   = *Pipes::DMEFT_results::Param["C63"];
      // float C64   = *Pipes::DMEFT_results::Param["C64"];
      // float C71   = *Pipes::DMEFT_results::Param["C71"]; 
      // float C72   = *Pipes::DMEFT_results::Param["C72"];
      // float C73   = *Pipes::DMEFT_results::Param["C73"];
      // float C74   = *Pipes::DMEFT_results::Param["C74"];
      // float mchi  = *Pipes::DMEFT_results::Param["mchi"]; 
      // float lambda= *Pipes::DMEFT_results::Param["Lambda"];

      // **--------------------------------------------------------------------------------------------//
      //** --------------------------------CMS---------------------------------------------------------//

      // Test the function to see if it compiles. 

      // const int CMS_SIZE = 22;
      #define CMS_SIZE 22
      double _srnums_CMS[CMS_SIZE];

      L_Acc_Eff_CS(_srnums_CMS, mchi,C61,C62,C63,C64,C71,C72,C73,C74,lambda,"CMS");
      
      // cout << "first _srnums call ..."<<endl;

      std::vector<SignalRegionData> cmsBinnedResults;
      
      // for (size_t ibin = 0; ibin < cms_bin_size; ++ibin) {
      //     std::stringstream ss; ss << "sr-" << ibin;
      //     cmsBinnedResults.push_back(SignalRegionData(ss.str(), OBSNUM[ibin], {_srnums[ibin],  0.}, {BKGNUM[ibin], BKGERR[ibin]}, _srnums[ibin]));
      // }
      
      for (size_t ibin = 0; ibin < cms_bin_size; ++ibin) {
        // Generate an 'sr-N' label 
        std::stringstream ss; ss << "sr-" << ibin;

        // Construct a SignalRegionData instance and add it to cmsBinnedResults
        SignalRegionData sr;
        sr.sr_label = ss.str();
        sr.n_observed = CMS_OBSNUM[ibin];
        sr.n_signal = _srnums_CMS[ibin];
        sr.n_signal_at_lumi = _srnums_CMS[ibin];  // We have already scaled the signals in _srnums_CMS to xsec * lumi
        sr.signal_sys = 0.;
        sr.n_background = CMS_BKGNUM[ibin];
        sr.background_sys = CMS_BKGERR[ibin];
        cmsBinnedResults.push_back(sr);
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

      static bool first_c = true;
      if (first_c){
        cout << " Defining CMS covariance matrix..."<<endl;
        for (int i = 0; i < 22; i++)
            m_BKGCOV.row(i) = Eigen::VectorXd::Map(&BKGCOV[i][0],BKGCOV[i].size()); 
        first_c = false;
      }

      // Save the results + covariance matrix in cmsData
      cmsData.srdata = cmsBinnedResults;
      cmsData.srcov = m_BKGCOV;

      // Save a copy of the results *without* the covariance matrix in cmsData_nocovar
      cmsData_nocovar.srdata = cmsBinnedResults;


      // **----------------------------------------------------------------------------------------------------//
      // **-------------------------------------ATLAS----------------------------------------------------------//

      // Now put the ATLAS data into an equivalent object
      // Andre to add the relevant lines
      

      // std::cout << "Making signal numbers" << std::endl; 

      // cout<<"Just b4 atlas srnums"<<endl;

      double _srnums_ATLAS[atlas_bin_size];
      L_Acc_Eff_CS(_srnums_ATLAS,mchi,C61,C62,C63,C64,C71,C72,C73,C74,lambda,"ATLAS"); 

      // cout << "Atlas srnums defined" <<endl;



      // cout << "After static atlas" <<endl;

      std::vector<SignalRegionData> atlasBinnedResults;

      for (size_t ibin = 0; ibin < atlas_bin_size; ++ibin) {

        // Generate an 'sr-N' label 
        std::stringstream ss; ss << "sr-" << ibin;
        // Construct a SignalRegionData instance and add it to atlasBinnedResults
        SignalRegionData sr;
        sr.sr_label = ss.str();
        sr.n_observed = ATLAS_OBSNUM[ibin];
        sr.n_signal = _srnums_ATLAS[ibin];
        sr.n_signal_at_lumi = _srnums_ATLAS[ibin];  // We have already scaled the signals in _srnums_ATLAS to xsec * lumi
        // cout << "Check output: "<< sr.sr_label<< "  " << _srnums_ATLAS[ibin] <<endl;
        sr.signal_sys = 0.;
        sr.n_background = ATLAS_BKGNUM[ibin];
        sr.background_sys = ATLAS_BKGERR[ibin];
        atlasBinnedResults.push_back(sr);
      }

      // Save the results in atlasData
      atlasData.srdata = atlasBinnedResults;


      // ******** Create total results ***********// 
      // //--------------------------------------//

      // Saving the addresses to the thread_local AnalysisData instances.
      result.push_back(&atlasData);
      result.push_back(&cmsData);
      result.push_back(&cmsData_nocovar);

      //Sleep time
      // std::this_thread::sleep_for(std::chrono::seconds(1));

      // auto finish_wall_clock = std::chrono::steady_clock::now();

      // Calculating total time taken by the program. 
      // cout << fixed << setprecision(8) << "Excecution time for DMEFT_results: " << ((finish_wall_clock - start_wall_clock) / std::chrono::nanoseconds(1))/(1E9) << '\n';
    };
         

    void calc_DMEFT_profiled_LHC_nuisance_params(map_str_dbl& result)
    {
      using namespace Pipes::calc_DMEFT_profiled_LHC_nuisance_params;

      // Steal the list of skipped analyses from the options from the "calc_combined_LHC_LogLike" function
      std::vector<str> default_skip_analyses;  // The default is empty lists of analyses to skip
      static const std::vector<str> skip_analyses = Pipes::calc_combined_LHC_LogLike::runOptions->getValueOrDef<std::vector<str> >(default_skip_analyses, "skip_analyses");
      
      // Clear previous result map
      result.clear();

      // Get Lambda
      const Spectrum& spec = *Dep::DMEFT_spectrum;
      float lambda = spec.get(Par::mass1, "Lambda");

      // 
      // Do profiling
      // 

      // Currently this is just a silly grid scan in 'a'. Will replace by a GLS minimizer
      float a = 1.0;
      float delta_a = 0.1;
      float a_bestfit = a;
      double total_loglike_bestfit = -50000;

      AnalysisLogLikes analoglikes;

      while (a < 100.0)
      {
        double total_loglike = 0.0;

        // Create a vector with temp AnalysisData instances by copying the original ones
        vector<AnalysisData> temp_adata_vec;
        for (AnalysisData* adata_ptr : *Dep::AllAnalysisNumbersUnmodified)
        {

          const string& analysis_name = adata_ptr->analysis_name;

          // If the analysis name is in skip_analyses, don't take it into account in this profiling
          if (std::find(skip_analyses.begin(), skip_analyses.end(), analysis_name) != skip_analyses.end())
          {
            cout << "DEBUG: Leaving out analysis " << analysis_name << " from profiling" << endl;
            continue;
          }

          // Make a copy of the AnalysisData instance that adata_ptr points to
          temp_adata_vec.push_back( AnalysisData(*adata_ptr) );
        }

        // Now loop over all the temp AnalysisData instances and calculate the total loglike for the current a-value
        for (AnalysisData& adata : temp_adata_vec)
        {
          const string& analysis_name = adata.analysis_name;

          signal_modifier_function(adata, lambda, a);
          analoglikes = calc_loglikes_for_analysis(adata, true, false, false);
          cout << "DEBUG: " << analysis_name << ": a, analoglikes.combination_loglike : " << a << ", " << analoglikes.combination_loglike << endl;

          total_loglike += analoglikes.combination_loglike;
        }

        // Save values if this is the current best-fit a-value
        if (total_loglike > total_loglike_bestfit)
        {
          total_loglike_bestfit = total_loglike;
          a_bestfit = a;
        }

        cout << "DEBUG: a, total_loglike: " << a << ", " << total_loglike << endl;

        // Next step
        a += delta_a;
      }

      cout << "DEBUG: " << endl;
      cout << "DEBUG: a_bestfit: " << a_bestfit << endl;
      cout << "DEBUG: total_loglike_bestfit: " << total_loglike_bestfit << endl;

      // Save the best-fit parameter value
      result["a"] = a_bestfit;
    }


    void DMEFT_results_profiled(AnalysisDataPointers& result)
    {
      using namespace Pipes::DMEFT_results_profiled;

      // Clear previous vectors, etc.
      result.clear();

      // Get the original AnalysisDataPointers that we will adjust
      result = *Dep::AllAnalysisNumbersUnmodified;

      // Get the best-fit nuisance parameter(s)
      map_str_dbl bestfit_nuisance_pars = *Dep::DMEFT_profiled_LHC_nuisance_params;
      float a_bestfit = bestfit_nuisance_pars.at("a");

      // Get Lambda
      const Spectrum& spec = *Dep::DMEFT_spectrum;
      float lambda = spec.get(Par::mass1, "Lambda");

      // Recalculate AnalysisData instances in "result", using the best-fit a-value
      for (AnalysisData* adata_ptr : result)
      {
        signal_modifier_function(*adata_ptr, lambda, a_bestfit);
      }
    }

   
    void InterpolatedMCInfo(MCLoopInfo& result)
    {
      // cout << "Have run the void..."<<endl;
      // This makes an MCLoopInfo object for satisfying the LHC
      // likelihood calculation dependency
      // ------------------------------- //

      // Event generation has been bypassed:
      // ------------------------------------------------------//
      result.event_gen_BYPASS = true;
      // ------------------------------------------------------//
      result.reset_flags();
      
    }
    
  } // namespace ColliderBit
    
} // namespace Gambit
