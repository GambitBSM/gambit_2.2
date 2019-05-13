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
#include "gambit/ColliderBit/lhef2heputils.hpp"
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
  // Declare a function that reads the total cross-section from the input file, but builds up the number of events from the event loop
  #define CAPABILITY CrossSection
    #define FUNCTION getYAMLxsec
    START_FUNCTION(xsec)
    #undef FUNCTION
  #undef CAPABILITY
#undef MODULE

namespace Gambit
{
  namespace ColliderBit
  {

    /// A nested function that reads in Les Houches Event files and converts them to HEPUtils::Event format
    void getLHEvent(HEPUtils::Event& result)
    {
      using namespace Pipes::getLHEvent;

      result.clear();

      // Get the pointer to the LHEF reader object
      LHEF::Reader& lhe = runOptions->getValue<LHEF::Reader&>("lhef_reader");

      // Don't do anything during special iterations
      cout << *Loop::iteration << endl;
      if (*Loop::iteration <= 0) return;

      // Attempt to read the next LHE event as a HEPUtils event. If there are no more events, wrap up the loop.
      #pragma omp critical
      {
        if (not lhe.readEvent()) Loop::wrapup();
        result = get_HEPUtils_event(lhe);
      }
    }

    /// A function that reads the total cross-section from the input file, but builds up the number of events from the event loop
    void getYAMLxsec(xsec& result)
    {
      using namespace Pipes::getYAMLxsec;

      // Don't bother if there are no analyses that will use this.
      if (Dep::RunMC->analyses.empty()) return;

      // Retrieve the total cross-section and cross-section error
      const static double xsec_fb = 1000 * runOptions->getValue<double>("xsec_pb");
      const static double xsec_fractional_uncert = runOptions->getValue<double>("xsec_fractional_uncert");

      // Reset the xsec objects on all threads
      if (*Loop::iteration == START_SUBPROCESS)
      {
        result.reset();
      }

      // If we are in the main event loop, count the event towards cross-section normalisation on this thread
      if (*Loop::iteration > 0)
      {
        result.log_event();
      }

      // Set the xsec and its error
      if (*Loop::iteration == COLLIDER_FINALIZE)
      {
        result.set_xsec(xsec_fb, xsec_fractional_uncert*xsec_fb);
        result.gather_num_events();
        cout << "xsec details: " << result.num_events() << " " << xsec_fb << " " << xsec_fractional_uncert << endl;
      }
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
    str lhef_filename = settings.getValue<str>("event_file");

    // Initialise logs
    logger().set_log_debug_messages(debug);
    initialise_standalone_logs("CBS_logs/");
    logger()<<"Running CBS"<<LogTags::info<<EOM;

    // Pass options to the main event loop
    YAML::Node CBS(infile["settings"]);
    CBS["analyses"] = analyses;
    CBS["min_nEvents"] = (long long)(1 000);
    CBS["max_nEvents"] = (long long)(1 000 000 000);
    operateLHCLoop.setOption<YAML::Node>("CBS", CBS);

    // Open the LHE file
    LHEF::Reader lhe(lhef_filename);

    // Pass (a pointer to) the LHEF reader object to the module function that uses it to return HEPUtils events
    getLHEvent.setOption<LHEF::Reader&>("lhef_reader", lhe);

    // Pass options to the cross-section function
    getYAMLxsec.setOption<double>("xsec_pb", settings.getValue<double>("xsec_pb"));
    getYAMLxsec.setOption<double>("xsec_fractional_uncert", settings.getValue<double>("xsec_fractional_uncert"));

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
    getATLASAnalysisContainer.resolveDependency(&getYAMLxsec);
    getCMSAnalysisContainer.resolveDependency(&getYAMLxsec);
    getIdentityAnalysisContainer.resolveDependency(&getYAMLxsec);
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
    getYAMLxsec.resolveLoopManager(&operateLHCLoop);
    runATLASAnalyses.resolveLoopManager(&operateLHCLoop);
    runCMSAnalyses.resolveLoopManager(&operateLHCLoop);
    runIdentityAnalyses.resolveLoopManager(&operateLHCLoop);
    vector<functor*> nested_functions = initVector<functor*>(&getLHEvent,
                                                             &getBuckFastATLAS,
                                                             &getBuckFastCMS,
                                                             &getBuckFastIdentity,
                                                             &getYAMLxsec,
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
    operateLHCLoop.reset_and_calculate();
    CollectAnalyses.reset_and_calculate();
    calc_LHC_LogLikes.reset_and_calculate();
    get_LHC_LogLike_per_analysis.reset_and_calculate();
    calc_combined_LHC_LogLike.reset_and_calculate();

    // Retrieve and print the predicted + observed counts and likelihoods for the individual SRs and analyses, as well as the total likelihood.
    /// @todo
    /// Get predicted counts (BG, signal, total), observed counts, likelihood for each SR
    /// Get selected SR information
    /// Get total likelihood for each analysis
    double loglike = calc_combined_LHC_LogLike(0);

    cout.precision(5);
    cout << endl;
    /// @todo
    /// print predicted counts (BG, signal, total), observed counts, likelihood for each SR
    /// print selected SR information
    /// print total likelihood for each analysis
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

