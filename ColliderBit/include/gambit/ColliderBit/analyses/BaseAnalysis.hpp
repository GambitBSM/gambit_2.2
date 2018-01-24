#pragma once
//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  The BaseAnalysis class.
///
/// @todo Move some inlines into .cpp files to minimise rebuilding?

#include "gambit/ColliderBit/ColliderBit_macros.hpp"
#include "gambit/ColliderBit/analyses/AnalysisData.hpp"

#include "gambit/ColliderBit/Utils.hpp"
#include "HEPUtils/MathUtils.h"
#include "HEPUtils/Event.h"

#include <string>
#include <sstream>
#include <vector>
#include <cmath>
#include <cfloat>
#include <limits>
#include <memory>
#include <iomanip>
#include <algorithm>

namespace Gambit {
  namespace ColliderBit {

    /// An abstract base class for collider analyses within ColliderBit.
    /// @note The templating makes forward declaration / #include decoupling hard :-(
    /// @todo How is this templating useful? Only if alternative Analysis implementations have the same interface...
    template <typename EventT>
    class BaseAnalysis {
    private:

      double _ntot, _xsec, _xsecerr, _luminosity;
      bool _xsec_is_set, _luminosity_is_set, _is_scaled;
      bool _needs_collection;
      AnalysisData _results;
      typedef EventT EventType;


    public:

      /// @name Construction, Destruction, and Recycling:
      //@{
      BaseAnalysis() : _ntot(0), _xsec(0), _xsecerr(0), _luminosity(0),
                       _xsec_is_set(false), _luminosity_is_set(false),
                       _is_scaled(false), _needs_collection(true) {  }
      virtual ~BaseAnalysis() { }
      /// Reset this instance for reuse, avoiding the need for "new" or "delete".
      virtual void clear() {
        _ntot = 0; _xsec = 0; _xsecerr = 0; _luminosity = 0; 
        _xsec_is_set = false; _luminosity_is_set = false;
        _is_scaled = false; _needs_collection = true;
        _results.clear();
      }
      //@}


      /// @name Event analysis, event number, and cross section functions:
      //@{
      /// Analyze the event (accessed by reference).
      void do_analysis(const EventT& e) { do_analysis(&e); }
      /// Analyze the event (accessed by pointer).
      void do_analysis(const EventT* e) { _needs_collection = true; analyze(e); }

      /// Return the total number of events seen so far.
      double num_events() const { return _ntot; }
      /// Return the cross-section (in pb).
      double xsec() const { return _xsec; }
      /// Return the cross-section error (in pb).
      double xsec_err() const { return _xsecerr; }
      /// Return the cross-section relative error.
      double xsec_relerr() const { return xsec() > 0 ? xsec_err()/xsec() : 0; }
      /// Return the cross-section per event seen (in pb).
      double xsec_per_event() const { return (xsec() >= 0 && num_events() > 0) ? xsec()/num_events() : 0; }
      /// Return the integrated luminosity (in inverse pb).
      double luminosity() const { return _luminosity; }
      /// Set the cross-section and its error (in pb).
      void set_xsec(double xs, double xserr) { _xsec_is_set = true; _xsec = xs; _xsecerr = xserr; }
      /// Set the integrated luminosity (in inverse pb).
      void set_luminosity(double lumi) { _luminosity_is_set = true; _luminosity = lumi; }

      /// Get the collection of SignalRegionData for likelihood computation.
      const AnalysisData& get_results()
      {
        if (_needs_collection)
        {
          collect_results();
          _needs_collection = false;
        }
        return _results;
      }

      /// An overload of get_results() with some additional consistency checks.
      const AnalysisData& get_results(std::string& warning)
      {
        warning = "";
        if (not _xsec_is_set)
          warning += "Cross section has not been set. ";
        if (not _luminosity_is_set)
          warning += "Luminosity has not been set. ";
        if (not _is_scaled)
          warning += "Results have not been scaled. ";
        if (_ntot < 1)
          warning += "No events have been analyzed. ";

        /// @todo We need to shift the 'analysis_name' property from class SignalRegionData 
        ///       to this class. Then we can add the class name to this error message.
        // warning = "Ooops! In analysis " + analysis_name + ": " + warning

        return get_results();
      }

      //@}


    protected:

      /// @name Protected collection functions
      //@{
      /// Analyze the event (accessed by pointer).
      /// @note Needs to be called from Derived::analyze().
      virtual void analyze(const EventT*) { _ntot += 1; }
      /// Add the given result to the internal results list.
      void add_result(const SignalRegionData& sr) { _results.add(sr); }
      /// Set the covariance matrix, expressing SR correlations
      void set_covariance(const Eigen::MatrixXd& srcov) { _results.srcov = srcov; }
      /// A convenience function for setting the SR covariance from a nested vector/initialiser list
      void set_covariance(const std::vector<std::vector<double>>& srcov) {
        Eigen::MatrixXd cov(srcov.size(), srcov.front().size());
        for (size_t i = 0; i < srcov.size(); ++i) {
          for (size_t j = 0; j < srcov.front().size(); ++j) {
            cov(i,j) = srcov[i][j];
          }
        }
        set_covariance(cov);
      }
      /// Gather together the info for likelihood calculation.
      virtual void collect_results() = 0;
      //@}


    public:

      /// @name (Re-)initialization functions
      //@{
      /// General init for any analysis of this type.
      virtual void init(const std::vector<std::string>&) {}
      /// General init for any collider of this type - no settings version.
      virtual void init() { }
      /// Scale by number of input events and xsec.
      virtual void scale(double factor=-1) {
        if (factor < 0) {
          factor = (num_events() == 0 ? 0 : (luminosity() * xsec()) / num_events());
          // cout << "DEBUG: " << luminosity() << " * " << xsec() << " / " << num_events() << " = " << factor << endl;
        }
        assert(factor >= 0);
        for (SignalRegionData& sr : _results) {
          sr.n_signal_at_lumi = factor * sr.n_signal;
          //cout << "DEBUG: " << factor << ", " << sr.n_signal << " -> " << sr.n_signal_at_lumi << endl;
        }
        _is_scaled = true;
      }
      //@}


      /// @name BaseAnalysis combination operations
      //@{
      /// An operator to do xsec-weighted combination of analysis runs.
      virtual void add(BaseAnalysis* other) {
        if (_results.empty()) collect_results();
        AnalysisData otherResults = other->get_results();
        /// @todo Access by name, including merging disjoint region sets?
        assert(otherResults.size() == _results.size());
        for (size_t i = 0; i < _results.size(); ++i) {
          _results[i].n_signal += otherResults[i].n_signal;
        }
        _ntot += other->num_events();
      }

      /// Add cross-sections and errors for two different process types.
      void add_xsec(double xs, double xserr) {
        if (xs > 0) {
          if (xsec() <= 0) {
            set_xsec(xs, xserr);
          } else {
            _xsec += xs;
            _xsecerr = HEPUtils::add_quad(xsec_err(), xserr);
          }
        }
      }
      /// Combine cross-sections and errors for the same process type, assuming uncorrelated errors.
      void improve_xsec(double xs, double xserr) {
        if (xs > 0) {
          if (xsec() <= 0) {
            set_xsec(xs, xserr);
          } else {
            /// @todo Probably shouldn't be combined with equal weight?!?
            _xsec = _xsec/2.0 + xs/2.0;
            _xsecerr = HEPUtils::add_quad(xsec_err(), xserr) / 2.0;
          }
        }
      }
      //@}

    };


    /// A BaseAnalysis template specialization for the standard event type.
    using HEPUtilsAnalysis = BaseAnalysis<HEPUtils::Event>;

  }
}
