//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  ColliderBit LHC signal and likelihood functions.
///
///  *********************************************
///
///  Authors (add name and date if you modify):
///
///  \author Abram Krislock
///          (a.m.b.krislock@fys.uio.no)
///
///  \author Aldo Saavedra
///
///  \author Andy Buckley
///
///  \author Chris Rogan
///          (crogan@cern.ch)
///  \date 2014 Aug
///  \date 2015 May
///
///  \author Pat Scott
///          (p.scott@imperial.ac.uk)
///  \date 2015 Jul
///  \date 2018 Jan
///  \date 2019 Jan
///
///  \author Anders Kvellestad
///          (anders.kvellestad@fys.uio.no)
///  \date   2017 March
///  \date   2018 Jan
///  \date   2018 May
///
///  *********************************************

#include <string>
#include <sstream>

#include "gambit/Elements/gambit_module_headers.hpp"
#include "gambit/ColliderBit/ColliderBit_rollcall.hpp"

#include "Eigen/Eigenvalues"
#include <gsl/gsl_sf_gamma.h>
#include "gambit/ColliderBit/multimin.h"

// #define COLLIDERBIT_DEBUG

namespace Gambit
{

  namespace ColliderBit
  {

    /// Loop over all analyses and fill a map of predicted counts
    void calc_LHC_signals(map_str_dbl& result)
    {
      using namespace Pipes::calc_LHC_signals;

      // Clear the result map
      result.clear();

      std::stringstream summary_line;
      summary_line << "LHC signals per SR: ";

      // Loop over analyses and collect the predicted events into the map
      for (size_t analysis = 0; analysis < Dep::AllAnalysisNumbers->size(); ++analysis)
      {
        // AnalysisData for this analysis
        const AnalysisData& adata = *(Dep::AllAnalysisNumbers->at(analysis));

        summary_line << adata.analysis_name << ": ";

        // Loop over the signal regions inside the analysis, and save the predicted number of events for each.
        for (size_t SR = 0; SR < adata.size(); ++SR)
        {
          // Save SR numbers and absolute uncertainties
          const SignalRegionData srData = adata[SR];
          const str key = adata.analysis_name + "__" + srData.sr_label + "__i" + std::to_string(SR) + "__signal";
          result[key] = srData.n_signal_at_lumi;
          const double abs_uncertainty_s_stat = (srData.n_signal == 0 ? 0 : sqrt(srData.n_signal) * (srData.n_signal_at_lumi/srData.n_signal));
          const double abs_uncertainty_s_sys = srData.signal_sys;
          const double combined_uncertainty = HEPUtils::add_quad(abs_uncertainty_s_stat, abs_uncertainty_s_sys);
          result[key + "_uncert"] = combined_uncertainty;

          summary_line << srData.sr_label + "__i" + std::to_string(SR) << ":" << srData.n_signal_at_lumi << "+-" << combined_uncertainty << ", ";
        }
      }
      logger() << LogTags::debug << summary_line.str() << EOM;
    }




    /// Loglike objective-function wrapper to provide the signature for GSL multimin
    ///
    /// @todo Use this to compute both the b and s+b likelihood, and compute
    /// the b likelihood only once, via s = 0
    void _gsl_calc_Analysis_MinusLogLike(const size_t n, const double* unit_nuisances_dbl,
                                         void* npred_nobs_cov, double* fval) {

      // Convert the array of doubles into an "Eigen view" of the nuisance params
      Eigen::Map<const Eigen::VectorXd> unit_nuisances(&unit_nuisances_dbl[0], n);

      // Convert the linearised array of doubles into "Eigen views" of the fixed params
      double *npred_nobs_cov_dbl = (double*) npred_nobs_cov;
      Eigen::Map<const Eigen::ArrayXd> n_pred_nominal(&npred_nobs_cov_dbl[0], n);
      Eigen::Map<const Eigen::ArrayXd> n_obs(&npred_nobs_cov_dbl[n], n);
      Eigen::Map<const Eigen::MatrixXd> cov(&npred_nobs_cov_dbl[2*n], n, n);

      // Call the real loglike function
      // return calc_Analysis_LogLikes(unit_nuisances, n_pred_nominal, n_obs, cov);

      // Rotate rate deltas into the SR basis and shift by SR mean rates
      const Eigen::VectorXd n_pred = n_pred_nominal + (cov*unit_nuisances).array();

      // Calculate each SR's Poisson likelihood and add to composite likelihood calculation
      double loglike_tot = 0;
      for (int j = 0; j < unit_nuisances.size(); ++j) {
        // First the multivariate Gaussian bit (j = nuisance)
        const double pnorm_j = -pow(unit_nuisances(j), 2)/2. - log(1/sqrt(2*M_PI));
        loglike_tot += pnorm_j;
        // Then the Poisson bit (j = SR)
        /// @note We've dropped the log(n_obs!) terms, since they're expensive and cancel in computing DLL
        const double lambda_j = std::max(n_pred(j), 1e-3); //< manually avoid <= 0 rates
        const double logfact_n_obs = 0; // gsl_sf_lngamma(n_obs(j) + 1); //< skipping log(n_obs!) computation
        const double loglike_j = n_obs(j)*log(lambda_j) - lambda_j - logfact_n_obs;
        loglike_tot += loglike_j;
      }

      // Output via argument (invert to return -LL for minimisation)
      *fval = -loglike_tot;
    }


    std::vector<double> _gsl_mkpackedarray(const Eigen::ArrayXd& n_pred,
                                           const Eigen::ArrayXd& n_obs,
                                           const Eigen::MatrixXd& cov) {
      const size_t nSR = n_obs.size();
      std::vector<double> fixeds(nSR + nSR + nSR*nSR, 0.0);
      for (size_t i = 0; i < nSR; ++i) {
        fixeds[0+i] = n_pred(i);
        fixeds[nSR+i] = n_obs(i);
        for (size_t j = 0; j < nSR; ++j)
          fixeds[2*nSR+i*nSR+j] = cov(i,j); ///< @todo Double-check ordering
      }
      return fixeds;
    }


    /// Return the best log likelihood
    /// @note Return value is missing the log(n_obs!) terms (n_SR of them) which cancel in LLR calculation
    double profile_loglike_cov(const Eigen::ArrayXd& n_pred,
                               const Eigen::ArrayXd& n_obs,
                               const Eigen::MatrixXd& cov) {

      // Number of signal regions
      const size_t nSR = n_obs.size();

      // Start with nuisances at nominal values
      /// @todo Pick a more informed starting position
      std::vector<double> nuisances(nSR, 0.0);

      // Optimiser parameters
      // Params: step1 size, tol, maxiter, epsabs, simplex maxsize, method, verbosity
      // Methods:
      //  0: Fletcher-Reeves conjugate gradient
      //  1: Polak-Ribiere conjugate gradient
      //  2: Vector Broyden-Fletcher-Goldfarb-Shanno method
      //  3: Steepest descent algorithm
      //  4: Nelder-Mead simplex
      //  5: Vector Broyden-Fletcher-Goldfarb-Shanno method ver. 2
      //  6: Simplex algorithm of Nelder and Mead ver. 2
      //  7: Simplex algorithm of Nelder and Mead: random initialization
      //
      /// @todo Tune / take from YAML: currently using 0.1 initial step, 0.01 convergence, and Simplex2
      struct multimin_params oparams = {.1,1e-2,100,1e-3,1e-5,6,0};

      // Convert the linearised array of doubles into "Eigen views" of the fixed params
      std::vector<double> fixeds = _gsl_mkpackedarray(n_pred, n_obs, cov);

      // Pass to the minimiser
      double minusbestll = 999;
      multimin(nSR, &nuisances[0], &minusbestll,
               nullptr, nullptr, nullptr,
               _gsl_calc_Analysis_MinusLogLike, nullptr, nullptr,
               &fixeds, oparams);

      return -minusbestll;
    }



    double marg_loglike_cov(const Eigen::ArrayXd& n_pred,
                            const Eigen::ArrayXd& n_obs,
                            const Eigen::ArrayXd& principlesd,
                            const Eigen::MatrixXd& cov) {

      // Number of signal regions
      const size_t nSR = n_obs.size();

      // Sample correlated SR rates from a rotated Gaussian defined by the covariance matrix and offset by the mean rates
      using namespace Pipes::calc_LHC_LogLikes;
      static const double CONVERGENCE_TOLERANCE_ABS = runOptions->getValueOrDef<double>(0.05, "covariance_marg_convthres_abs");
      static const double CONVERGENCE_TOLERANCE_REL = runOptions->getValueOrDef<double>(0.05, "covariance_marg_convthres_rel");
      static const size_t NSAMPLE_INPUT = runOptions->getValueOrDef<size_t>(100000, "covariance_nsamples_start");

      // Dynamic convergence control & test variables
      size_t nsample = NSAMPLE_INPUT;
      bool first_iteration = true;
      double diff_abs = 9999;
      double diff_rel = 1;

      // Likelihood variables (note use of long double to guard against blow-up of L as opposed to log(L1/L0))
      long double ana_like_prev = 1;
      long double ana_like = 1;
      long double lsum_prev = 0;

      // Sampler for unit-normal nuisances
      std::normal_distribution<double> unitnormdbn(0,1);

      // Log factorial of observed number of events.
      // Currently use the ln(Gamma(x)) function gsl_sf_lngamma from GSL. (Need continuous function.)
      // We may want to switch to using Stirling's approximation: ln(n!) ~ n*ln(n) - n
      Eigen::ArrayXd logfact_n_obs(nSR);
      for (size_t j = 0; j < nSR; ++j)
        logfact_n_obs(j) = gsl_sf_lngamma(n_obs(j) + 1);

      // Check absolute difference between independent estimates
      /// @todo Should also implement a check of relative difference
      while ((diff_abs > CONVERGENCE_TOLERANCE_ABS && diff_rel > CONVERGENCE_TOLERANCE_REL) || 1.0/sqrt(nsample) > CONVERGENCE_TOLERANCE_ABS)
      {
        long double lsum = 0;

        /// @note How to correct negative rates? Discard (scales badly), set to
        /// epsilon (= discontinuous & unphysical pdf), transform to log-space
        /// (distorts the pdf quite badly), or something else (skew term)?
        /// We're using the "set to epsilon" version for now.
        /// Ben: I would vote for 'discard'. It can't be that inefficient, surely?
        /// Andy: For a lot of signal regions, the probability of none having a negative sample is Prod_SR p(non-negative)_SR... which *can* get bad.

        #pragma omp parallel
        {

          ////////////////////
          // Start one-point likelihood calculation
          /// @todo Chop out into a separate function? Just the bit in the for loop, to keep OMP happy?
          ////////////////////

          // Sample correlated SR rates from a rotated Gaussian defined by the covariance matrix and offset by the mean rates
          double lsum_private  = 0;
          #pragma omp for nowait
          for (size_t i = 0; i < nsample; ++i) {

            Eigen::VectorXd norm_sample(nSR);
            for (size_t j = 0; j < nSR; ++j)
              norm_sample(j) = principlesd(j) * unitnormdbn(Random::rng());

            // Rotate rate deltas into the SR basis and shift by SR mean rates
            const Eigen::VectorXd n_pred_sample  = n_pred + (cov*norm_sample).array();

            // Calculate Poisson likelihood and add to composite likelihood calculation
            double combined_loglike = 0;
            for (size_t j = 0; j < nSR; ++j) {
              const double lambda_j = std::max(n_pred_sample(j), 1e-3); //< manually avoid <= 0 rates
              const double loglike_j  = n_obs(j)*log(lambda_j) - lambda_j - logfact_n_obs(j);
              combined_loglike  += loglike_j;
            }
            // Add combined likelihood to running sums (to later calculate averages)
            lsum_private += exp(combined_loglike);
          }

          ////////////////////
          /// ^^^ end one-point likelihood calculation
          ////////////////////

          #pragma omp critical
          {
            lsum  += lsum_private;
          }

        } // End omp parallel

        // Compare convergence to previous independent batch
        if (first_iteration)  // The first round must be generated twice
        {
          lsum_prev = lsum;
          first_iteration = false;
        }
        else
        {
          ana_like_prev = lsum_prev / (double)nsample;
          ana_like = lsum / (double)nsample;
          diff_abs = fabs(ana_like_prev - ana_like);
          diff_rel = diff_abs/ana_like;

          // Update variables
          lsum_prev += lsum;  // Aggregate result. This doubles the effective batch size for lsum_prev.
          nsample *=2;  // This ensures that the next batch for lsum is as big as the current batch size for lsum_prev, so they can be compared directly.
        }

        #ifdef COLLIDERBIT_DEBUG
        cout << debug_prefix()
             << "diff_rel: " << diff_rel << endl
             << "   diff_abs: " << diff_abs << endl
             << "   logl: " << log(ana_like) << endl;
        cout << debug_prefix() << "nsample for the next iteration is: " << nsample << endl;
        cout << debug_prefix() << endl;
        #endif
      }
      // End convergence while-loop

      // Combine the independent estimates ana_like and ana_like_prev.
      // Use equal weights since the estimates are based on equal batch sizes.
      ana_like = 0.5*(ana_like + ana_like_prev);
      const double ana_margll = log(ana_like);
      #ifdef COLLIDERBIT_DEBUG
      cout << debug_prefix() << "Combined estimate: ana_loglike: " << ana_margll << "   (based on 2*nsample=" << 2*nsample << " samples)" << endl;
      #endif

      return ana_margll;
    }




    /// Loop over all analyses and fill a map of AnalysisLogLikes objects
    void calc_LHC_LogLikes(map_str_AnalysisLogLikes& result)
    {
      using namespace Pipes::calc_LHC_LogLikes;

      // Use marginalisation rather than profiling (probably less stable)?
      static const bool USE_MARG = runOptions->getValueOrDef<bool>(false, "use_marginalising");

      // Use covariance matrices if available?
      static const bool USE_COVAR = runOptions->getValueOrDef<bool>(true, "use_covariances");

      // Clear the result map
      result.clear();

      // Main loop over all analyses to compute DLL = LL_sb - LL_b
      for (size_t analysis = 0; analysis < Dep::AllAnalysisNumbers->size(); ++analysis)
      {


        // AnalysisData for this analysis
        const AnalysisData& adata = *(Dep::AllAnalysisNumbers->at(analysis));

        /// If no events have been generated (xsec veto) or too many events have failed,
        /// short-circut the loop and return delta log-likelihood = 0 for every SR in
        /// each analysis.
        /// @todo This must be made more sophisticated once we add analyses that
        ///       don't rely on event generation.
        if (not Dep::RunMC->event_generation_began || Dep::RunMC->exceeded_maxFailedEvents)
        {
          // If this is an analysis with covariance info, only add a single 0-entry in the map
          if (USE_COVAR && adata.srcov.rows() > 0)
          {
            result[adata.analysis_name].combination_sr_label = "none";
            result[adata.analysis_name].combination_loglike = 0.0;
            continue;
          }
          // If this is an analysis without covariance info, add 0-entries for all SRs plus
          // one for the combined LogLike
          else
          {
            for (size_t SR = 0; SR < adata.size(); ++SR)
            {
              result[adata.analysis_name].sr_indices[adata[SR].sr_label] = SR;
              result[adata.analysis_name].sr_loglikes[adata[SR].sr_label] = 0.0;
              continue;
            }
            result[adata.analysis_name].combination_sr_label = "none";
            result[adata.analysis_name].combination_loglike = 0.0;
            continue;
          }

        }


        #ifdef COLLIDERBIT_DEBUG
        std::streamsize stream_precision = cout.precision();  // get current precision
        cout.precision(2);  // set precision
        cout << debug_prefix() << "calc_LHC_LogLikes: " << "Will print content of " << adata.analysis_name << " signal regions:" << endl;
        for (size_t SR = 0; SR < adata.size(); ++SR)
        {
          const SignalRegionData& srData = adata[SR];
          cout << std::fixed << debug_prefix()
                                 << "calc_LHC_LogLikes: " << adata.analysis_name
                                 << ", " << srData.sr_label
                                 << ",  n_b = " << srData.n_background << " +/- " << srData.background_sys
                                 << ",  n_obs = " << srData.n_observed
                                 << ",  excess = " << srData.n_observed - srData.n_background << " +/- " << srData.background_sys
                                 << ",  n_s = " << srData.n_signal_at_lumi
                                 << ",  (excess-n_s) = " << (srData.n_observed-srData.n_background) - srData.n_signal_at_lumi << " +/- " << srData.background_sys
                                 << ",  n_s_MC = " << srData.n_signal
                                 << endl;
        }
        cout.precision(stream_precision); // restore previous precision
        #endif


        // Loop over the signal regions inside the analysis, and work out the total (delta) log likelihood for this analysis
        /// @todo Unify the treatment of best-only and correlated SR treatments as far as possible
        /// @todo Come up with a good treatment of zero and negative predictions
        if (USE_COVAR && adata.srcov.rows() > 0)
        {
          /// If (simplified) SR-correlation info is available, so use the
          /// covariance matrix to construct composite marginalised likelihood
          /// Despite initial thoughts, we can't just do independent LL
          /// calculations in a rotated basis, but have to sample from the
          /// covariance matrix.
          ///
          /// @note This means we can't use the nulike LL functions, which
          /// operate in 1D only.  Also, log-normal sampling in the diagonal
          /// basis is not helpful, since the rotation will re-generate negative
          /// rates.
          ///
          /// @todo How about Gaussian sampling in the log(rate) space? Would
          /// protect against negatives in any SR. Requires care with the
          /// explicit transformation of widths.
          ///
          /// @todo Support skewness correction to the pdf.

          #ifdef COLLIDERBIT_DEBUG
          cout << debug_prefix() << "calc_LHC_LogLikes: Analysis " << analysis << " has a covariance matrix: computing composite loglike." << endl;
          #endif


          // Shortcut: if all SRs have 0 signal prediction, we know the Delta LogLike is 0.
          bool all_zero_signal = true;
          size_t nSR = adata.size();
          for (size_t SR = 0; SR < nSR; ++SR)
          {
            if (adata[SR].n_signal != 0)
            {
              all_zero_signal = false;
              break;
            }
          }
          if (all_zero_signal)
          {
            // Store result
            result[adata.analysis_name].combination_sr_label = "all";
            result[adata.analysis_name].combination_sr_index = -1;
            result[adata.analysis_name].combination_loglike = 0.0;

            #ifdef COLLIDERBIT_DEBUG
            cout << debug_prefix() << "calc_LHC_LogLikes: " << adata.analysis_name << "_LogLike : " << 0.0 << " (No signal predicted. Skipped covariance calculation.)" <<endl;
            #endif

            // Continue to next analysis
            continue;
          }

          // Construct vectors of SR numbers
          Eigen::ArrayXd n_obs(adata.size()), logfact_n_obs(adata.size()), n_pred_b(adata.size()), n_pred_sb(adata.size()), abs_unc_s(adata.size());
          for (size_t SR = 0; SR < adata.size(); ++SR)
          {
            const SignalRegionData srData = adata[SR];

            // Actual observed number of events
            n_obs(SR) = srData.n_observed;

            // Log factorial of observed number of events.
            // Currently use the ln(Gamma(x)) function gsl_sf_lngamma from GSL. (Need continuous function.)
            // We may want to switch to using Stirling's approximation: ln(n!) ~ n*ln(n) - n
            logfact_n_obs(SR) = gsl_sf_lngamma(n_obs(SR) + 1.);

            // A contribution to the predicted number of events that is not known exactly
            n_pred_b(SR) = srData.n_background;
            n_pred_sb(SR) = srData.n_signal_at_lumi + srData.n_background;

            // Absolute errors for n_predicted_uncertain_*
            const double abs_uncertainty_s_stat = (srData.n_signal == 0 ? 0 : sqrt(srData.n_signal) * (srData.n_signal_at_lumi/srData.n_signal));
            const double abs_uncertainty_s_sys = srData.signal_sys;
            abs_unc_s(SR) = HEPUtils::add_quad(abs_uncertainty_s_stat, abs_uncertainty_s_sys);
          }

          // Diagonalise the background-only covariance matrix, extracting the rotation matrix
          /// @todo No need to recompute the background-only covariance decomposition for every point!
          /// Ben: Actually don't need to recompute the background-only marginalisation at all. It
          ///      is always the same, so can just do it once at the start of the scan.
          const Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eig_b(adata.srcov);
          const Eigen::ArrayXd Eb = eig_b.eigenvalues();
          const Eigen::ArrayXd sqrtEb = Eb.sqrt();
          const Eigen::MatrixXd Vb = eig_b.eigenvectors();
          //const Eigen::MatrixXd Vbinv = Vb.inverse();

          // Construct and diagonalise the s+b covariance matrix, adding the diagonal signal uncertainties in quadrature
          /// @todo Is this the best way, or should we just sample the s numbers independently and then be able to completely cache the cov matrix diagonalisation?
          const Eigen::MatrixXd srcov_s = abs_unc_s.array().square().matrix().asDiagonal();
          const Eigen::MatrixXd srcov_sb = adata.srcov + srcov_s;
          const Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eig_sb(srcov_sb);
          const Eigen::ArrayXd Esb = eig_sb.eigenvalues();
          const Eigen::ArrayXd sqrtEsb = Esb.sqrt();
          const Eigen::MatrixXd Vsb = eig_sb.eigenvectors();
          //const Eigen::MatrixXd Vsbinv = Vsb.inverse();

          // Analysis-level DLL result
          double ana_dll = NAN;

          if (!USE_MARG) {

            // Pass background to the profiler
            /// @todo Only compute this once per run
            const double llpart_b = profile_loglike_cov(n_pred_b, n_obs, Vb);

            // Pass signal+background to the profiler
            const double llpart_sb = profile_loglike_cov(n_pred_sb, n_obs, Vsb);

            // Compute the DLL wrt background-only
            ana_dll = llpart_sb - llpart_b;

          } else {

            // Pass background to the marginaliser
            /// @todo Only compute this once per run
            const double ll_b = marg_loglike_cov(n_pred_b, n_obs, sqrtEb, Vb);

            // Pass signal+background to the marginaliser
            const double ll_sb = marg_loglike_cov(n_pred_sb, n_obs, sqrtEsb, Vsb);

            // Compute the DLL wrt background-only
            ana_dll = ll_sb - ll_b;


            ////////////////////
            /// start likelihood marginalisation
            /// @todo Split this whole chunk off into a lnlike-style utility function?
            ////////////////////

            // // Sample correlated SR rates from a rotated Gaussian defined by the covariance matrix and offset by the mean rates
            // static const double CONVERGENCE_TOLERANCE_ABS = runOptions->getValueOrDef<double>(0.05, "covariance_marg_convthres_abs");
            // static const double CONVERGENCE_TOLERANCE_REL = runOptions->getValueOrDef<double>(0.05, "covariance_marg_convthres_rel");
            // static const size_t NSAMPLE_INPUT = runOptions->getValueOrDef<size_t>(100000, "covariance_nsamples_start");
            // size_t nsample = NSAMPLE_INPUT;

            // // Dynamic convergence control & test variables
            // bool first_iteration = true;
            // double diff_abs = 9999;
            // double diff_rel = 1;

            // // Likelihood variables (note use of long double to guard against blow-up of L as opposed to log(L1/L0))
            // long double ana_like_b_prev = 1;
            // long double ana_like_sb_prev = 1;
            // long double ana_like_b = 1;
            // long double ana_like_sb = 1;
            // long double lsum_b_prev = 0;
            // long double lsum_sb_prev = 0;

            // std::normal_distribution<double> unitnormdbn(0,1);

            // // Check absolute difference between independent estimates
            // /// @todo Should also implement a check of relative difference
            // while ((diff_abs > CONVERGENCE_TOLERANCE_ABS && diff_rel > CONVERGENCE_TOLERANCE_REL) || 1.0/sqrt(nsample) > CONVERGENCE_TOLERANCE_ABS)
            // {
            //   long double lsum_b = 0;
            //   long double lsum_sb = 0;

            //   /// @note How to correct negative rates? Discard (scales badly), set to
            //   /// epsilon (= discontinuous & unphysical pdf), transform to log-space
            //   /// (distorts the pdf quite badly), or something else (skew term)?
            //   /// We're using the "set to epsilon" version for now.
            //   /// Ben: I would vote for 'discard'. It can't be that inefficient, surely?
            //   /// Andy: For a lot of signal regions, the probability of none having a negative sample is Prod_SR p_SR(non-negative)... which *can* get bad.

            //   #pragma omp parallel
            //   {

            //     ////////////////////
            //     /// start one-point likelihood calculation
            //     ////////////////////

            //     double lsum_b_private  = 0;
            //     double lsum_sb_private = 0;

            //     // Sample correlated SR rates from a rotated Gaussian defined by the covariance matrix and offset by the mean rates
            //     #pragma omp for nowait
            //     for (size_t i = 0; i < nsample; ++i) {

            //       Eigen::VectorXd norm_sample_b(adata.size()), norm_sample_sb(adata.size());
            //       for (size_t j = 0; j < adata.size(); ++j) {
            //         norm_sample_b(j) = sqrtEb(j) * unitnormdbn(Random::rng());
            //         norm_sample_sb(j) = sqrtEsb(j) * unitnormdbn(Random::rng());
            //       }

            //       // Rotate rate deltas into the SR basis and shift by SR mean rates
            //       const Eigen::VectorXd n_pred_b_sample  = n_pred_b + (Vb*norm_sample_b).array();
            //       const Eigen::VectorXd n_pred_sb_sample = n_pred_sb + (Vsb*norm_sample_sb).array();

            //       // Calculate Poisson likelihood and add to composite likelihood calculation
            //       double combined_loglike_b = 0;
            //       double combined_loglike_sb = 0;
            //       for (size_t j = 0; j < adata.size(); ++j) {
            //         const double lambda_b_j = std::max(n_pred_b_sample(j), 1e-3); //< manually avoid <= 0 rates
            //         const double lambda_sb_j = std::max(n_pred_sb_sample(j), 1e-3); //< manually avoid <= 0 rates
            //         const double loglike_b_j  = n_obs(j)*log(lambda_b_j) - lambda_b_j - logfact_n_obs(j);
            //         const double loglike_sb_j = n_obs(j)*log(lambda_sb_j) - lambda_sb_j - logfact_n_obs(j);
            //         combined_loglike_b  += loglike_b_j;
            //         combined_loglike_sb += loglike_sb_j;
            //       }
            //       // Add combined likelihood to running sums (to later calculate averages)
            //       lsum_b_private  += exp(combined_loglike_b);
            //       lsum_sb_private += exp(combined_loglike_sb);
            //     }

            //     ////////////////////
            //     /// ^^^ end one-point likelihood calculation
            //     ////////////////////

            //     #pragma omp critical
            //     {
            //       lsum_b  += lsum_b_private;
            //       lsum_sb += lsum_sb_private;
            //     }

            //   } // End omp parallel

            //   // Compare convergence to previous independent batch
            //   if (first_iteration)  // The first round must be generated twice
            //   {
            //     lsum_b_prev = lsum_b;
            //     lsum_sb_prev = lsum_sb;
            //     first_iteration = false;
            //   }
            //   else
            //   {
            //     ana_like_b_prev = lsum_b_prev / (double)nsample;
            //     ana_like_sb_prev = lsum_sb_prev / (double)nsample;
            //     ana_like_b = lsum_b / (double)nsample;
            //     ana_like_sb = lsum_sb / (double)nsample;
            //     //
            //     const double diff_abs_b = fabs(ana_like_b_prev - ana_like_b);
            //     const double diff_abs_sb = fabs(ana_like_sb_prev - ana_like_sb);
            //     const double diff_rel_b = diff_abs_b/ana_like_b;
            //     const double diff_rel_sb = diff_abs_sb/ana_like_sb;
            //     //
            //     diff_rel = std::max(diff_rel_b, diff_rel_sb);  // Relative convergence check
            //     diff_abs = std::max(diff_abs_b, diff_abs_sb);  // Absolute convergence check

            //     // Update variables
            //     lsum_b_prev += lsum_b;  // Aggregate result. This doubles the effective batch size for lsum_prev.
            //     lsum_sb_prev += lsum_sb;  // Aggregate result. This doubles the effective batch size for lsum_prev.
            //     nsample *=2;  // This ensures that the next batch for lsum is as big as the current batch size for lsum_prev, so they can be compared directly.
            //   }

            //   #ifdef COLLIDERBIT_DEBUG
            //   cout << debug_prefix()
            //        << "diff_rel: " << diff_rel << endl
            //        <<  "   diff_abs: " << diff_abs << endl
            //        << "   ana_llr_prev: " << log(ana_like_sb_prev/ana_like_b_prev) << endl
            //        << "   ana_dll: " << log(ana_like_sb/ana_like_b) << endl
            //        << "   logl_sb: " << log(ana_like_sb) << endl
            //        << "   logl_b: " << log(ana_like_b) << endl;
            //   cout << debug_prefix() << "nsample for the next iteration is: " << nsample << endl;
            //   cout << debug_prefix() << endl;
            //   #endif
            // }  // End while loop

            // // Combine the independent estimates ana_like and ana_like_prev.
            // // Use equal weights since the estimates are based on equal batch sizes.
            // ana_like_b = 0.5*(ana_like_b + ana_like_b_prev);
            // ana_like_sb = 0.5*(ana_like_sb + ana_like_sb_prev);

            // // Compute LLR from mean s+b and b likelihoods
            // ana_dll = log(ana_like_sb) - log(ana_like_b);
            // #ifdef COLLIDERBIT_DEBUG
            // cout << debug_prefix() << "Combined estimate: ana_dll: " << ana_dll << "   (based on 2*nsample=" << 2*nsample << " samples)" << endl;
            // #endif

            ///////////////////////
            /// ^^^ end likelihood marginalisation
            ///////////////////////

          } // end ana_dll calculation, by either marginalisation or profiling


          // Check for problems with the result
          if (Utils::isnan(ana_dll))
          {
            std::stringstream msg;
            msg << "Computation of composite loglike for analysis " << adata.analysis_name << " returned NaN. Will now print the signal region data for this analysis:" << endl;
            for (size_t SR = 0; SR < adata.size(); ++SR)
            {
              const SignalRegionData& srData = adata[SR];
              msg << srData.sr_label
                  << ",  n_background = " << srData.n_background
                  << ",  background_sys = " << srData.background_sys
                  << ",  n_observed = " << srData.n_observed
                  << ",  n_signal_at_lumi = " << srData.n_signal_at_lumi
                  << ",  n_signal = " << srData.n_signal
                  << ",  signal_sys = " << srData.signal_sys
                  << endl;
            }
            invalid_point().raise(msg.str());
          }

          // Store result
          result[adata.analysis_name].combination_sr_label = "all";
          result[adata.analysis_name].combination_sr_index = -1;
          result[adata.analysis_name].combination_loglike = ana_dll;

          #ifdef COLLIDERBIT_DEBUG
          cout << debug_prefix() << "calc_LHC_LogLikes: " << adata.analysis_name << "_LogLike : " << ana_dll << endl;
          #endif

        }

        else
        {

          /// @todo Huge duplication here... reuse existing functions for trivial vectors & cov={sigma2}
          /// @todo Nulike is marginalising rather than profiling over the single best-exp SR: make consistent

          // No SR-correlation info, or user chose not to use it.
          // Then we either take the result from the SR *expected* to be most constraining
          // under the s=0 assumption (default), or naively combine the loglikes for
          // all SRs (if combine_SRs_without_covariances=true).
          #ifdef COLLIDERBIT_DEBUG
          cout << debug_prefix() << "calc_LHC_LogLikes: Analysis " << analysis << " has no covariance matrix: computing single best-expected loglike." << endl;
          #endif

          double bestexp_dll_exp = 0, bestexp_dll_obs = 0;
          str bestexp_sr_label;
          int bestexp_sr_index;
          double nocovar_srsum_dll_obs = 0;

          for (size_t SR = 0; SR < adata.size(); ++SR)
          {
            const SignalRegionData &srData = adata[SR];

            // Actual observed number of events
            const int n_obs = (int) round(srData.n_observed);

            // A contribution to the predicted number of events that is known exactly
            // (e.g. from data-driven background estimate)
            const double n_predicted_exact = 0;

            // A contribution to the predicted number of events that is not known exactly
            const double n_predicted_uncertain_s = srData.n_signal_at_lumi;
            const double n_predicted_uncertain_b = srData.n_background;
            const double n_predicted_uncertain_sb = n_predicted_uncertain_s + n_predicted_uncertain_b;

            // Absolute errors for n_predicted_uncertain_*
            const double abs_uncertainty_s_stat = (srData.n_signal == 0 ? 0 : sqrt(srData.n_signal) * (srData.n_signal_at_lumi/srData.n_signal));
            const double abs_uncertainty_s_sys = srData.signal_sys;
            const double abs_uncertainty_b = srData.background_sys;
            const double abs_uncertainty_sb = HEPUtils::add_quad(abs_uncertainty_s_stat, abs_uncertainty_s_sys, abs_uncertainty_b);

            // Relative errors for n_predicted_uncertain_*
            const double frac_uncertainty_b = abs_uncertainty_b / n_predicted_uncertain_b;
            const double frac_uncertainty_sb = abs_uncertainty_sb / n_predicted_uncertain_sb;

            // Predicted total background, as an integer for use in Poisson functions
            const int n_predicted_total_b_int = (int) round(n_predicted_exact + n_predicted_uncertain_b);

            // Marginalise over systematic uncertainties on mean rates
            // Use a log-normal/Gaussia distribution for the nuisance parameter, as requested
            auto marginaliser = (*BEgroup::lnlike_marg_poisson == "lnlike_marg_poisson_lognormal_error")
              ? BEreq::lnlike_marg_poisson_lognormal_error : BEreq::lnlike_marg_poisson_gaussian_error;
            const double llb_exp =  marginaliser(n_predicted_total_b_int, n_predicted_exact, n_predicted_uncertain_b, frac_uncertainty_b);
            const double llsb_exp = marginaliser(n_predicted_total_b_int, n_predicted_exact, n_predicted_uncertain_sb, frac_uncertainty_sb);
            const double llb_obs =  marginaliser(n_obs, n_predicted_exact, n_predicted_uncertain_b, frac_uncertainty_b);
            const double llsb_obs = marginaliser(n_obs, n_predicted_exact, n_predicted_uncertain_sb, frac_uncertainty_sb);

            // Calculate the expected dll and set the bestexp values for exp and obs dll if this one is the best so far
            const double dll_exp = llb_exp - llsb_exp; //< note positive dll convention -> more exclusion here
            #ifdef COLLIDERBIT_DEBUG
            cout << debug_prefix() << adata.analysis_name << ", " << srData.sr_label << ",  llsb_exp-llb_exp = " << llsb_exp-llb_exp << ",  llsb_obs-llb_obs= " << llsb_obs - llb_obs << endl;
            #endif
            if (dll_exp > bestexp_dll_exp || SR == 0)
            {
              bestexp_dll_exp = dll_exp;
              bestexp_dll_obs = llb_obs - llsb_obs;
              bestexp_sr_label = srData.sr_label;
              bestexp_sr_index = SR;
              // #ifdef COLLIDERBIT_DEBUG
              // cout << debug_prefix() << "Setting bestexp_sr_label to: " << bestexp_sr_label << ", LogL_exp = " << -bestexp_dll_exp << ", LogL_obs = " << -bestexp_dll_obs << endl;
              // #endif
            }

            // Store "observed LogLike" result for this SR
            result[adata.analysis_name].sr_indices[srData.sr_label] = SR;
            result[adata.analysis_name].sr_loglikes[srData.sr_label] = llsb_obs - llb_obs;

            // Add loglike to the no-correlations loglike sum over SRs
            nocovar_srsum_dll_obs += llsb_obs - llb_obs;
          }

          // Check for problem
          if (Utils::isnan(bestexp_dll_obs))
          {
            std::stringstream msg;
            msg << "Computation of single-SR loglike for analysis " << adata.analysis_name << " returned NaN, from signal region: " << bestexp_sr_label << endl;
            msg << "Will now print the signal region data for this analysis:" << endl;
            for (size_t SR = 0; SR < adata.size(); ++SR)
            {
              const SignalRegionData& srData = adata[SR];
              msg << srData.sr_label
                  << ",  n_background = " << srData.n_background
                  << ",  background_sys = " << srData.background_sys
                  << ",  n_observed = " << srData.n_observed
                  << ",  n_signal_at_lumi = " << srData.n_signal_at_lumi
                  << ",  n_signal = " << srData.n_signal
                  << ",  signal_sys = " << srData.signal_sys
                  << endl;
            }
            invalid_point().raise(msg.str());
          }

          // Set this analysis' total obs dLL to that from the best-expected SR (with conversion to more negative dll = more exclusion convention)
          // result[adata.analysis_name] = -bestexp_dll_obs;
          result[adata.analysis_name].combination_sr_label = bestexp_sr_label;
          result[adata.analysis_name].combination_sr_index = bestexp_sr_index;
          result[adata.analysis_name].combination_loglike = -bestexp_dll_obs;

          // Should we use the naive sum of SR loglikes (without correlations), instead of the best-expected SR?
          static const bool combine_nocovar_SRs = runOptions->getValueOrDef<bool>(false, "combine_SRs_without_covariances");
          if (combine_nocovar_SRs)
          {
            result[adata.analysis_name].combination_loglike = nocovar_srsum_dll_obs;
          }

          #ifdef COLLIDERBIT_DEBUG
          cout << debug_prefix() << "calc_LHC_LogLikes: " << adata.analysis_name << "_" << bestexp_sr_label << "_LogLike : " << -bestexp_dll_obs << endl;
          #endif
        }

      }

    }


    /// Extract the combined log likelihood for each analysis
    void get_LHC_LogLike_per_analysis(map_str_dbl& result)
    {
      using namespace Pipes::get_LHC_LogLike_per_analysis;

      std::stringstream summary_line;
      summary_line << "LHC loglikes per analysis: ";

      for (const std::pair<str,AnalysisLogLikes>& pair : *Dep::LHC_LogLikes)
      {
        const str& analysis_name = pair.first;
        const AnalysisLogLikes& analysis_loglikes = pair.second;

        result[analysis_name] = analysis_loglikes.combination_loglike;

        summary_line << analysis_name << ":" << analysis_loglikes.combination_loglike << ", ";
      }
      logger() << LogTags::debug << summary_line.str() << EOM;
    }


    /// Extract the log likelihood for each SR
    void get_LHC_LogLike_per_SR(map_str_dbl& result)
    {
      using namespace Pipes::get_LHC_LogLike_per_SR;

      std::stringstream summary_line;
      summary_line << "LHC loglikes per SR: ";

      for (const std::pair<str,AnalysisLogLikes>& pair_i : *Dep::LHC_LogLikes)
      {
        const str& analysis_name = pair_i.first;
        const AnalysisLogLikes& analysis_loglikes = pair_i.second;

        summary_line << analysis_name << ": ";

        for (const std::pair<str,double>& pair_j : analysis_loglikes.sr_loglikes)
        {
          const str& sr_label = pair_j.first;
          const double& sr_loglike = pair_j.second;
          const int sr_index = analysis_loglikes.sr_indices.at(sr_label);

          const str key = analysis_name + "__" + sr_label + "__i" + std::to_string(sr_index) + "__LogLike";
          result[key] = sr_loglike;

          summary_line << sr_label + "__i" + std::to_string(sr_index) << ":" << sr_loglike << ", ";
        }

        result[analysis_name + "__combined_LogLike"] = analysis_loglikes.combination_loglike;

        summary_line << "combined_LogLike:" << analysis_loglikes.combination_loglike << ", ";
      }
      logger() << LogTags::debug << summary_line.str() << EOM;
    }


    /// Extract the labels for the SRs used in the analysis loglikes
    void get_LHC_LogLike_SR_labels(map_str_str& result)
    {
      using namespace Pipes::get_LHC_LogLike_per_SR;
      for (const std::pair<str,AnalysisLogLikes>& pair_i : *Dep::LHC_LogLikes)
      {
        const str& analysis_name = pair_i.first;
        const AnalysisLogLikes& analysis_loglikes = pair_i.second;

        result[analysis_name] = analysis_loglikes.combination_sr_label;
      }
    }


    /// Extract the indices for the SRs used in the analysis loglikes
    /// @todo Switch result type to map_str_int once we have implemented a printer for this type
    void get_LHC_LogLike_SR_indices(map_str_dbl& result)
    {
      using namespace Pipes::get_LHC_LogLike_per_SR;

      std::stringstream summary_line;
      summary_line << "LHC loglike SR indices: ";

      // Loop over analyses
      for (const std::pair<str,AnalysisLogLikes>& pair_i : *Dep::LHC_LogLikes)
      {
        const str& analysis_name = pair_i.first;
        const AnalysisLogLikes& analysis_loglikes = pair_i.second;

        result[analysis_name] = (double) analysis_loglikes.combination_sr_index;

        summary_line << analysis_name << ":" << analysis_loglikes.combination_sr_index << ", ";
      }
      logger() << LogTags::debug << summary_line.str() << EOM;
    }


    /// Compute the total likelihood combining all analyses
    void calc_combined_LHC_LogLike(double& result)
    {
      using namespace Pipes::calc_combined_LHC_LogLike;
      result = 0.0;

      // Read analysis names from the yaml file
      std::vector<str> default_skip_analyses;  // The default is empty lists of analyses to skip
      static const std::vector<str> skip_analyses = runOptions->getValueOrDef<std::vector<str> >(default_skip_analyses, "skip_analyses");

      // If too many events have failed, do the conservative thing and return delta log-likelihood = 0
      if (Dep::RunMC->exceeded_maxFailedEvents)
      {
        #ifdef COLLIDERBIT_DEBUG
          cout << debug_prefix() << "calc_combined_LHC_LogLike: Too many failed events. Will be conservative and return a delta log-likelihood of 0." << endl;
        #endif
        return;
      }

      // Loop over analyses and calculate the total observed dLL
      for (auto const& analysis_loglike_pair : *Dep::LHC_LogLike_per_analysis)
      {
        const str& analysis_name = analysis_loglike_pair.first;
        const double& analysis_loglike = analysis_loglike_pair.second;

        // If the analysis name is in skip_analyses, don't add its loglike to the total loglike.
        if (std::find(skip_analyses.begin(), skip_analyses.end(), analysis_name) != skip_analyses.end())
        {
          #ifdef COLLIDERBIT_DEBUG
            cout.precision(5);
            cout << debug_prefix() << "calc_combined_LHC_LogLike: Leaving out analysis " << analysis_name << " with LogL = " << analysis_loglike << endl;
          #endif
          continue;
        }

        // Add analysis loglike.
        // If using capped likelihood for each individual analysis, set analysis_loglike = min(analysis_loglike,0)
        static const bool use_cap_loglike_individual = runOptions->getValueOrDef<bool>(false, "cap_loglike_individual_analyses");
        if (use_cap_loglike_individual)
        {
          result += std::min(analysis_loglike, 0.0);
        }
        else
        {
          result += analysis_loglike;
        }

        #ifdef COLLIDERBIT_DEBUG
          cout.precision(5);
          cout << debug_prefix() << "calc_combined_LHC_LogLike: Analysis " << analysis_name << " contributes with a LogL = " << analysis_loglike << endl;
        #endif
      }

      #ifdef COLLIDERBIT_DEBUG
        cout << debug_prefix() << "calc_combined_LHC_LogLike: LHC_Combined_LogLike = " << result << endl;
      #endif

      // If using a "global" capped likelihood, set result = min(result,0)
      static const bool use_cap_loglike = runOptions->getValueOrDef<bool>(false, "cap_loglike");
      if (use_cap_loglike)
      {
        result = std::min(result, 0.0);
      }

      std::stringstream summary_line;
      summary_line << "LHC combined loglike:" << result;
      logger() << LogTags::debug << summary_line.str() << EOM;
    }

  }

}
