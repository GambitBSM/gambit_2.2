//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  Container for EFT parameterisation of WIMP
///  annihilations to SM particles
///
///  *********************************************
///
///  Authors (add name and date if you modify):
///
///  \author Ben Farmer
///          (b.farmer@imperial.ac.uk)
///  \date 2019 Sep
///
///  *********************************************

#ifndef __wimp_annihilation_hpp__
#define __wimp_annihilation_hpp__

#include <string>
#include <map>

namespace Gambit
{
    /// Contain for generic parameterisation of WIMP annihilation to various two-body final states,
    /// with <sigma v> expanded as a simple power series in v^2
    class WIMP_annihilation
    {
        public:
            WIMP_annihilation();
            double A(const std::string& channel) const;
            double B(const std::string& channel) const;

            void setA(const std::string& channel, double val);
            void setB(const std::string& channel, double val);
        private:
            std::map<std::string,double> a;
            std::map<std::string,double> b;
    };
}

#endif //__wimp_annihilation_hpp__
