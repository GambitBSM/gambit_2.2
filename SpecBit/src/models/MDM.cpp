//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  Functions of module SpecBit
///
///  SpecBit module functions related to the
///  MDM model.
///
///  *********************************************
///
///  Authors (add name and date if you modify):
///
///  \author James McKay
///    \date 2018 March
///
///  *********************************************

#include <string>
#include <sstream>

#include "gambit/Elements/gambit_module_headers.hpp"
#include "gambit/Elements/spectrum.hpp"

#include "gambit/SpecBit/SpecBit_rollcall.hpp"
//#include "gambit/SpecBit/SpecBit_helpers.hpp"
#include "gambit/SpecBit/RegisteredSpectra.hpp"

// Switch for debug mode
//#define SPECBIT_DEBUG

namespace Gambit
{

  namespace SpecBit
  {

    template <class T>
    void fill_MDM_input(T& input, const std::map<str, safe_ptr<const double> >& Param,SMInputs sminputs)
    {
      double mH = *Param.at("mH");
      double mChi = *Param.at("mChi");

      //double QEFT  = mChi;
			/*
      if (QEFT < sminputs.mT)
      {
				QEFT = sminputs.mT;
			}
      */
      double QEFT  = sminputs.mT;

      input.HiggsIN = 0.5*pow(mH,2);

      input.YcIN = 0.5*mChi;

      input.QEWSB = QEFT;  // scale where EWSB conditions are applied
      input.Qin = QEFT; //scale;  // highest scale at which model is run to

    }

    /* TODO: Not ready yet
    bool check_perturb_MDM(const Spectrum& spec,double scale,int pts)
    {
      std::unique_ptr<SubSpectrum> MDM = spec.clone_HE();
      double step = log10(scale) / pts;
      double runto;

      //const double ul = std::sqrt(4.0 * pi); // Maximum value for perturbative couplings, same perturbativity bound that FlexibleSUSY uses
      double ul = 4.0 * pi;
      for (int i=0;i<pts;i++)
      {
        runto = pow(10,step*float(i+1.0)); // scale to run spectrum to
        if (runto<100){runto=100.0;}// avoid running to low scales

        try
	      {
	        MDM -> RunToScale(runto);
	      }
	      catch (const Error& error)
	      {
	        return false;
	      };




        static const SpectrumContents::MDM contents;
        static const std::vector<SpectrumParameter> required_parameters = contents.all_parameters_with_tag(Par::dimensionless);

        for(std::vector<SpectrumParameter>::const_iterator it = required_parameters.begin();
              it != required_parameters.end(); ++it)
        {
          const Par::Tags        tag   = it->tag();
          const std::string      name  = it->name();
          const std::vector<int> shape = it->shape();
          std::ostringstream label;
          label << name <<" "<< Par::toString.at(tag);

          if (name == "lambda_h"){ul =  2*pi;}
          else {ul = 4.0 * pi;}

          if(shape.size()==1 and shape[0]==1)
          {
            if (abs(MDM->get(tag,name))>ul)
            {
							return false;
						}
          }
          else if(shape.size()==1 and shape[0]>1)
          {
            for(int k = 1; k<=shape[0]; ++k)
            {
              if (abs(MDM->get(tag,name,k))>ul) return false;
            }
          }
          else if(shape.size()==2)
          {
            for(int k = 1; k<=shape[0]; ++k)
            {
              for(int j = 1; j<=shape[0]; ++j)
              {
                if (abs(MDM->get(tag,name,k,j))>ul) return false;
              }
            }
          }

        }
      }

      return true;

    }*/

    void get_MDM_spectrum(Spectrum& result)
    {
      namespace myPipe = Pipes::get_MDM_spectrum;
      const SMInputs& sminputs = *myPipe::Dep::SMINPUTS;
      const Options& runOptions=*myPipe::runOptions;

      // Set up the input structure
      SpectrumInputs inputs(sminputs, SpectrumContents::MDM(),  myPipe::Param, myPipe::runOptions);

      // TODO: This is handled by the backend
      //fill_MDM_input(input,myPipe::Param,sminputs);

      // Get the spectrum from the Backend
      myPipe::BEreq::FS_MDM_Spectrum(result, inputs);

      // Retrieve any mass cuts
      result.check_mass_cuts(*myPipe::runOptions);

      int check_perturb_pts = runOptions.getValueOrDef<double>(10,"check_perturb_pts");
      double do_check_perturb = runOptions.getValueOrDef<bool>(false,"check_perturb");
      double check_perturb_scale = runOptions.getValueOrDef<double>(1.22e19,"check_high_scale");
      double input_scale_tolerance = runOptions.getValueOrDef<double>(1e-3,"input_scale_tolerance");

      // Check that the Higgs MSbar parameter has been provided at the scale of the MDM MSbar multiplet mass
      double Qin = *myPipe::Param.at("Qin");
      if (abs((Qin - *myPipe::Param.at("mChi"))/Qin) > input_scale_tolerance) SpecBit_error().raise(LOCAL_INFO, "SM_Higgs_running::Qin must equal MDM::mChi");

      if (do_check_perturb)
      {
        /* TODO: Not ready yet
        if (!check_perturb_MDM(result,check_perturb_scale,check_perturb_pts))
        {
          // invalidate point as spectrum not perturbative up to scale
          std::ostringstream msg;
          msg << "Spectrum not perturbative up to scale = " << check_perturb_scale <<  std::endl;
          #ifdef SPECBIT_DEBUG
            cout << "Spectrum not perturbative up to scale = " << check_perturb_scale <<  endl;
          #endif
          invalid_point().raise(msg.str());
        }*/
      }


    }


    void find_non_perturb_scale_MDM(double &result)
    {

      namespace myPipe = Pipes::find_non_perturb_scale_MDM;

      const Spectrum& fullspectrum = *myPipe::Dep::MDM_spectrum;

      // bound x by (a,b)
      // do all this is log space please

      double ms = *myPipe::Param.at("mS");

      double a = log10(ms);

      if (a > 20.0)
      {
        std::ostringstream msg;
        msg << "Scalar mass larger than 10^20 GeV " << std::endl;
        invalid_point().raise(msg.str());
      }

      double b = 20.0;
      double x = 0.5 * ( b + ms );

      while (abs(a-b)>1e-10)
      {
        //cout<< "\r" << "(a,b) = " << a << "  " << b << endl;
        //std::cout << std::flush;
        x=0.5*(b-a)+a;

        /* TODO: Not ready yet
        if (!check_perturb_MDM(fullspectrum,pow(10,x),3))
        {
          b=x;
        }
        else
        {
          a=x;
        }*/
      }
      result = pow(10,0.5*(a+b));
    }


    /// Print MDM spectrum out. Stripped down copy from MSSM version with variable names changed
    void fill_map_from_MDMspectrum(std::map<std::string,double>&, const Spectrum&);

    void get_MDM_spectrum_as_map (std::map<std::string,double>& specmap)
    {
      namespace myPipe = Pipes::get_MDM_spectrum_as_map;
      const Spectrum& mdmspec(*myPipe::Dep::MDM_spectrum);
      fill_map_from_MDMspectrum(specmap, mdmspec);
    }

    void fill_map_from_MDMspectrum(std::map<std::string,double>& specmap, const Spectrum& mdmspec)
    {
      /// Add everything... use spectrum contents routines to automate task
      static const SpectrumContents::MDM contents;
      static const std::vector<SpectrumContents::Parameter> required_parameters = contents.all_parameters();

      for(std::vector<SpectrumContents::Parameter>::const_iterator it = required_parameters.begin();
           it != required_parameters.end(); ++it)
      {
         const Par::Tags        tag   = it->tag();
         const std::string      name  = it->name();
         const std::vector<int> shape = it->shape();

         /// Verification routine should have taken care of invalid shapes etc, so won't check for that here.

         // Check scalar case
         if(shape.size()==1 and shape[0]==1)
         {
           std::ostringstream label;
           label << name <<" "<< Par::toString.at(tag);
           specmap[label.str()] = mdmspec.get(tag,name);
         }
         // Check vector case
         else if(shape.size()==1 and shape[0]>1)
         {
           for(int i = 1; i<=shape[0]; ++i) {
             std::ostringstream label;
             label << name <<"_"<<i<<" "<< Par::toString.at(tag);
             specmap[label.str()] = mdmspec.get(tag,name,i);
           }
         }
         // Check matrix case
         else if(shape.size()==2)
         {
           for(int i = 1; i<=shape[0]; ++i) {
             for(int j = 1; j<=shape[0]; ++j) {
               std::ostringstream label;
               label << name <<"_("<<i<<","<<j<<") "<<Par::toString.at(tag);
               specmap[label.str()] = mdmspec.get(tag,name,i,j);
             }
           }
         }
         // Deal with all other cases
         else
         {
           // ERROR
           std::ostringstream errmsg;
           errmsg << "Error, invalid parameter received while converting MDMspectrum to map of strings! This should no be possible if the spectrum content verification routines were working correctly; they must be buggy, please report this.";
           errmsg << "Problematic parameter was: "<< tag <<", " << name << ", shape="<< shape;
           utils_error().forced_throw(LOCAL_INFO,errmsg.str());
         }
      }

    }

  } // end namespace SpecBit
} // end namespace Gambit

