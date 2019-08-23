//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  ColliderBit (production) cross-section class.
///
///  *********************************************
///
///  Authors (add name and date if you modify):
///
///  \author Pat Scott
///          (p.scott@imperial.ac.uk)
///  \date 2019 Feb
///
///  *********************************************

#include <map>
#include <string>

#pragma once

namespace Gambit
{

  namespace ColliderBit
  {

    /// A class for holding cross-section info within ColliderBit.
    class xsec
    {

      public:

        xsec();
        virtual ~xsec() { }

        /// Reset this instance for reuse.
        void reset();

        /// Tell the xsec object that there has been a new event.
        void log_event();

        /// Return the total number of events seen so far.
        long long num_events() const;

        /// Return the full cross-section (in pb).
        double operator()() const;

        /// Return the cross-section error (in pb).
        double xsec_err() const;

        /// Return the cross-section relative error.
        double xsec_relerr() const;

        /// Return the cross-section per event seen (in pb).
        double xsec_per_event() const;

        /// Set the total number of events seen so far.
        void set_num_events(long long);

        /// Set the cross-section and its error (in pb).
        void set_xsec(double, double);

        /// Average cross-sections and combine errors.
        void average_xsec(double, double, long long);

        /// Collect xsec predictions from other threads and do a weighted combination.
        void gather_xsecs();

        /// Collect total events seen on all threads.
        void gather_num_events();

        /// Get content as map <string,double> map (for easy printing).
        std::map<std::string, double> get_content_as_map() const;

      private:

        long long _ntot;
        double _xsec;
        double _xsecerr;

        /// A map with pointers to all instances of this class. The key is the OMP thread number.
        static std::map<int, const xsec*> instances_map;

    };

  }
}
