//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  Class for holding model parameters. Defines the
///  basic container and get/set functions for
///  the model parameters. Originally adapted from
///  SUfit version (the connection is almost gone 
///  now though.)
///
///  *********************************************
///
///  Authors (add name and date if you modify):
///   
///  \author Johan Lundberg (SUfit version; model_types.hpp)
///  \date 2011 Jul - Aug
///
///  \author Ben Farmer
///  \date 2013 May, Jun 
///  \date 2014 Aug:  Went through and cut out some 
///                   stuff that seemed unused: removed
///                   the old ModelParametersBase class
///                   that didn't do much; removed the
///                   virtualisation since nothing derives
///                   from this class at the moment; 
///                   incorporated the proper error 
///                   handling system; put class 
///                   definitions into a seperate source file.
///
///  \author Pat Scott  
///          (patscott@physics.mcgill.ca)
///  \date 2013 Oct
///  \date 2014 Jan, Nov
///
///  *********************************************


#ifndef __model_parameters_hpp__
#define __model_parameters_hpp__

#include <map>
#include <vector>
#include <iostream>
#include <sstream>

#include "gambit/Utils/standalone_error_handlers.hpp"

namespace Gambit {

  // Model parameter map type; used by all models
  typedef std::map<std::string, double> parameterMap;

  class ModelParameters
  {

    protected:

      /// Checks if this model container holds a parameter match the supplied name
      void assert_contains(std::string) const;

    public:

      /// Default constructor
      ModelParameters();

      /// Constructor using vector of strings
      ModelParameters(const std::vector<std::string>&);
    
      /// Constructor using array of char arrays
      ModelParameters(const char**);
   
      /// Get value of named parameter 
      double getValue(std::string const & inkey) const;
    
      /// Get values of all parameters
      std::map<std::string, double> getValues() const;
   
      /// Get pointer to map of all parameters 
      ///TODO(Ben) I was going to delete this, but it seems to be used in some macros somewhere. 
      ///          I haven't checked what for, but probably whatever the macros want this pointer
      ///          for can be replaced by a new member function, to maintain better encapsulation.
      const std::map<std::string, double>* getValuesPtr() const; 

      /// Get number of parameters stored in this object
      int getNumberOfPars() const;

      /// Get parameter value using bracket operator
      const double & operator[](std::string const & inkey) const;

      /// Set single parameter value
      void setValue(std::string const &inkey,double const&value);
  
      /// Set many parameter values using a map
      void setValues(std::map<std::string,double> const &params_map);
    
      /// Get parameter keys (names), probably for external iteration
      std::vector<std::string> getKeys() const;

      /// Dump parameter names and values to stdout (should be for debugging only)
      void print() const;

      /// Dump parameter names and values to stream (again for debugging only I think)
      friend std::ostream &operator<<(std::ostream &strm, const ModelParameters &me)
      {
        strm << "ModelParameters: Printing: "<<std::endl;
        for (std::map<std::string,double>::const_iterator it=me._values.begin();it!=me._values.end();it++)
        {
          strm << "parameter: " << it->first << " value: "<<it->second ;
        }
        return strm;
      }

      /// Define a parameter with name, value (i.e. add to internal map). Value is initialised to zero
      void _definePar(const std::string &newkey);

      /// Define many new parameters at once via a vector of names
      void _definePars(const std::vector<std::string> &v);

      /// Define many new parameters at once via an array of char arrays
      void _definePars(const char** array);
  
    private:

      /// Internal map representation of parameters and their values
      std::map<std::string,double> _values;

  };

}
#endif //ifdef __model_parameters_hpp__