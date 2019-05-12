//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///
///  ColliderBit Solo: an event-based LHC recast
///  tool using the GAMBIT ColliderBit module.
///
///  *********************************************
///
///  Authors (add name and date if you modify):
///
///  \author Pat Scott
///  (p.scott@imperial.ac.uk)
///  \date May 2019
///
///  *********************************************

using namespace std;

#include <algorithm>
#include <dirent.h>

#include "LHEF.h"

#include "gambit/Elements/standalone_module.hpp"
#include "gambit/ColliderBit/ColliderBit_rollcall.hpp"
#include "gambit/Utils/cats.hpp"
#include "gambit/Utils/util_functions.hpp"

#define NULIKE_VERSION "1.0.7"
#define NULIKE_SAFE_VERSION 1_0_7

using namespace ColliderBit::Functown;
using namespace BackendIniBit::Functown;
using namespace CAT(Backends::nulike_,NULIKE_SAFE_VERSION)::Functown;

#define MODULE ColliderBit
  // Declare a nested function that reads in Les Houches Event files and converts them to HEPUtils::Event format
  #define CAPABILITY HardScatteringEvent
    #define FUNCTION getLHEvent
    START_FUNCTION(HEPUtils::Event)
    NEEDS_MANAGER(RunMC, MCLoopInfo)
    #undef FUNCTION
  #undef CAPABILITY
  // Declare a non-nested function that reads the total cross-section from the input file
  #define CAPABILITY CrossSection
    #define FUNCTION getXsecFromYaml
    START_FUNCTION(xsec)
    #undef FUNCTION
  #undef CAPABILITY
#undef MODULE

/// A nested function that reads in Les Houches Event files and converts them to HEPUtils::Event format
namespace Gambit
{
  namespace ColliderBit
  {
    void getLHEvent(HEPUtils::Event& result)
    {
      using namespace Pipes::getLHEvent;
      result.clear();
      // Retrieve the total number of LHEF files to be read
      static const long long n_lhefs = runOptions->getValue<long long>("n_lhefs");
      // Don't do anything during special iterations
      cout << *Loop::iteration << endl;
      if (*Loop::iteration <= 0) return;
      // Skip this loop iteration if we are out of LHEF files
      if (*Loop::iteration > n_lhefs) Loop::cycle();
      // Access the file corresponding to this iteration
      static const std::vector<str> filenames = runOptions->getValue<std::vector<str>>("lhef_filenames");
      str filename = filenames[*Loop::iteration-1];
      cout << "files: " << endl;
      for (auto x : filenames) {cout << "  " << x << endl;}
      cout << "Would have read file " << *Loop::iteration << " of " << n_lhefs << " on thread " << omp_get_thread_num() << ": " << filename << endl;
      // Convert the file into a HEPUtils::Event
      //result = read_LHEF(filename);
    }
  }
}

/// A non-nested function that reads the total cross-section from the input file
namespace Gambit
{
  namespace ColliderBit
  {
    void getXsecFromYaml(xsec& result)
    {
      using namespace Pipes::getXsecFromYaml;
      // Retrieve the total number of LHEF files to be read, total cross-section and cross-section error
      const long long n_lhefs = runOptions->getValue<long long>("n_lhefs");
      const double xsec_pb = runOptions->getValue<double>("xsec_pb");
      const double xsec_fractional_uncert = runOptions->getValue<double>("xsec_fractional_uncert");
      // Add them to the xsec object
      result.set_num_events(n_lhefs);
      result.set_xsec(xsec_pb, xsec_fractional_uncert*xsec_pb);
      cout << "xsec details: " << n_lhefs << " " << xsec_pb << " " << xsec_fractional_uncert << endl;
    }
  }
}


/// CBS main program
int main(int argc, char* argv[])
{

  try
  {

    // Check the number of command line arguments
    if (argc < 2)
    {
      // Tell the user how to run the program and exit
      cerr << endl << "Usage: " << argv[0] << " <your CBS yaml file>" << endl << endl;
      return 1;
    }

    // Make sure that nulike is present.
    if (not Backends::backendInfo().works[str("nulike")+NULIKE_VERSION]) backend_error().raise(LOCAL_INFO, str("nulike ")+NULIKE_VERSION+" is missing!");

    // Print the banner (if you could call it that)
    cout << endl;
    cout << "==================================" << endl;
    cout << "||                              ||" << endl;
    cout << "||    CBS: ColliderBit Solo     ||" << endl;
    cout << "||  GAMBIT Collider Workgroup   ||" << endl;
    cout << "||                              ||" << endl;
    cout << "==================================" << endl;
    cout << endl;

    // Read input file name
    const string filename_in = argv[1];

    // Read the settings in the input file
    YAML::Node infile;
    std::vector<str> analyses;
    Options settings;
    try
    {
      // Load up the file
      infile = YAML::LoadFile(filename_in);
      // Retrieve the analyses
      if (infile["analyses"]) analyses = infile["analyses"].as<std::vector<str>>();
      else throw std::runtime_error("Analyses list not found in "+filename_in+".  Quitting...");
      // Retrieve the other settings
      if (infile["settings"]) settings = Options(infile["settings"]);
      else throw std::runtime_error("Settings section not found in "+filename_in+".  Quitting...");
    }
    catch (YAML::Exception &e)
    {
      throw std::runtime_error("YAML error in "+filename_in+".\n(yaml-cpp error: "+std::string(e.what())+" )");
    }

    // Translate relevant settings into appropriate variables
    bool debug = settings.getValueOrDef<bool>(false, "debug");
    bool use_lnpiln = settings.getValueOrDef<bool>(false, "use_lognormal_distribution_for_1d_systematic");
    str lhef_dir = settings.getValue<str>("event_file_dir");
    DIR* dp = opendir(lhef_dir.c_str());
    if (dp == NULL) throw std::runtime_error("LHE file directory "+lhef_dir+" not found.  Quitting...");
    std::vector<str> lhef_filenames = Gambit::Utils::ls_dir(lhef_dir);
    long long n_lhefs = lhef_filenames.size();

    // Initialise logs
    logger().set_log_debug_messages(debug);
    initialise_standalone_logs("CBS_logs/");
    logger()<<"Running CBS"<<LogTags::info<<EOM;

    // Pass options to the main event loop
    YAML::Node CBS(infile["settings"]);
    CBS["analyses"] = analyses;
    CBS["min_nEvents"] = std::min((long long)(1000), n_lhefs);
    CBS["max_nEvents"] = n_lhefs;
    operateLHCLoop.setOption<YAML::Node>("CBS", CBS);

    // Pass options to the LHEF reader
    getLHEvent.setOption<long long>("n_lhefs", n_lhefs);
    getLHEvent.setOption<std::vector<str>>("lhef_filenames", lhef_filenames);

    // Pass options to the cross-section function
    getXsecFromYaml.setOption<long long>("n_lhefs", n_lhefs);
    getXsecFromYaml.setOption<double>("xsec_pb", settings.getValue<double>("xsec_pb"));
    getXsecFromYaml.setOption<double>("xsec_fractional_uncert", settings.getValue<double>("xsec_fractional_uncert"));

    // Pass options to the likelihood function
    calc_LHC_LogLikes.setOption<int>("covariance_nsamples_start", settings.getValue<int>("covariance_nsamples_start"));
    calc_LHC_LogLikes.setOption<double>("covariance_marg_convthres_abs", settings.getValue<double>("covariance_marg_convthres_abs"));
    calc_LHC_LogLikes.setOption<double>("covariance_marg_convthres_rel", settings.getValue<double>("covariance_marg_convthres_rel"));

    // Resolve ColliderBit dependencies and backend requirements
    calc_combined_LHC_LogLike.resolveDependency(&get_LHC_LogLike_per_analysis);
    calc_combined_LHC_LogLike.resolveDependency(&operateLHCLoop);
    get_LHC_LogLike_per_analysis.resolveDependency(&calc_LHC_LogLikes);
    calc_LHC_LogLikes.resolveDependency(&CollectAnalyses);
    calc_LHC_LogLikes.resolveDependency(&operateLHCLoop);
    calc_LHC_LogLikes.resolveBackendReq(use_lnpiln ? &nulike_lnpiln : &nulike_lnpin);
    CollectAnalyses.resolveDependency(&runATLASAnalyses);
    CollectAnalyses.resolveDependency(&runCMSAnalyses);
    CollectAnalyses.resolveDependency(&runIdentityAnalyses);
    runATLASAnalyses.resolveDependency(&getATLASAnalysisContainer);
    runATLASAnalyses.resolveDependency(&smearEventATLAS);
    runCMSAnalyses.resolveDependency(&getCMSAnalysisContainer);
    runCMSAnalyses.resolveDependency(&smearEventCMS);
    runIdentityAnalyses.resolveDependency(&getIdentityAnalysisContainer);
    runIdentityAnalyses.resolveDependency(&copyEvent);
    getATLASAnalysisContainer.resolveDependency(&getXsecFromYaml);
    getCMSAnalysisContainer.resolveDependency(&getXsecFromYaml);
    getIdentityAnalysisContainer.resolveDependency(&getXsecFromYaml);
    smearEventATLAS.resolveDependency(&getBuckFastATLAS);
    smearEventATLAS.resolveDependency(&getLHEvent);
    smearEventCMS.resolveDependency(&getBuckFastCMS);
    smearEventCMS.resolveDependency(&getLHEvent);
    copyEvent.resolveDependency(&getBuckFastIdentity);
    copyEvent.resolveDependency(&getLHEvent);

    // Resolve loop manager for ColliderBit event loop
    getLHEvent.resolveLoopManager(&operateLHCLoop);
    getBuckFastATLAS.resolveLoopManager(&operateLHCLoop);
    getBuckFastCMS.resolveLoopManager(&operateLHCLoop);
    getBuckFastIdentity.resolveLoopManager(&operateLHCLoop);
    getATLASAnalysisContainer.resolveLoopManager(&operateLHCLoop);
    getCMSAnalysisContainer.resolveLoopManager(&operateLHCLoop);
    getIdentityAnalysisContainer.resolveLoopManager(&operateLHCLoop);
    smearEventATLAS.resolveLoopManager(&operateLHCLoop);
    smearEventCMS.resolveLoopManager(&operateLHCLoop);
    copyEvent.resolveLoopManager(&operateLHCLoop);
    runATLASAnalyses.resolveLoopManager(&operateLHCLoop);
    runCMSAnalyses.resolveLoopManager(&operateLHCLoop);
    runIdentityAnalyses.resolveLoopManager(&operateLHCLoop);
    vector<functor*> nested_functions = initVector<functor*>(&getLHEvent,
                                                             &getBuckFastATLAS,
                                                             &getBuckFastCMS,
                                                             &getBuckFastIdentity,
                                                             &getATLASAnalysisContainer,
                                                             &getCMSAnalysisContainer,
                                                             &getIdentityAnalysisContainer,
                                                             &smearEventATLAS,
                                                             &smearEventCMS,
                                                             &copyEvent,
                                                             &runATLASAnalyses,
                                                             &runCMSAnalyses,
                                                             &runIdentityAnalyses);
    operateLHCLoop.setNestedList(nested_functions);

    // Call the initialisation function for nulike
    nulike_1_0_7_init.reset_and_calculate();

    // Run the detector sim and selected analyses on all the events read in.
    getXsecFromYaml.reset_and_calculate();
    operateLHCLoop.reset_and_calculate();
    CollectAnalyses.reset_and_calculate();
    calc_LHC_LogLikes.reset_and_calculate();
    get_LHC_LogLike_per_analysis.reset_and_calculate();
    calc_combined_LHC_LogLike.reset_and_calculate();

    // Retrieve and print the predicted + observed counts and likelihoods for the individual SRs and analyses, as well as the total likelihood.
    // get predicted counts (BG, signal, total), observed counts, likelihood for each SR
    // get selected SR information
    // get total likelihood for each analysis
    double loglike = calc_combined_LHC_LogLike(0);

    cout.precision(5);
    cout << endl;
    // print predicted counts (BG, signal, total), observed counts, likelihood for each SR
    // print selected SR information
    // print total likelihood for each analysis
    cout << scientific << "Total combined ATLAS+CMS log-likelihood: " << loglike << endl;
    cout << endl;

    // No more to see here folks, go home.
    return 0;
  }

  catch (std::exception& e)
  {
    cerr << "CBS has exited with fatal exception: " << e.what() << endl;
  }

  // Finished, but an exception was raised.
  return 1;

}

