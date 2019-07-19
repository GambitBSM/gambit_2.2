//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  Variadic utilty functions which work with
///  YAML objects.
///
///  Separated from variadic_functions.hpp to
///  avoid having the YAML headers included
///  everywhere.
///
///  *********************************************
///
///  Authors (add name and date if you modify):
///
///  \author Gregory Martinez
///          (gregory.david.martinez@gmail.com)
///  \date Feb 2014
//
///  \author Christoph Weniger
///          <c.weniger@uva.nl>
///  \date Dec 2014
///
///  \author Ben Farmer
///          <b.farmer@imperial.ac.uk>
///  \date Jan 2016, Jul 2019
///
///  *********************************************

#ifndef YAML_VARIADIC_FUNCTIONS_HPP
#define YAML_VARIADIC_FUNCTIONS_HPP

#include <string>

#include "yaml-cpp/yaml.h"
#include "gambit/Utils/standalone_error_handlers.hpp"

namespace Gambit
{       
        //////////////////////////////////////
        //Variadic Node functions
        //////////////////////////////////////
        
        inline const YAML::Node getVariadicNode(const YAML::Node &node)
        {
                return node;
        }
        
        inline const YAML::Node getVariadicNode(const YAML::Node &node, std::string key) 
        {
                if(node.Type()!=YAML::NodeType::Map)
                {
                     std::ostringstream os;
                     os<<"Error accessing YAML key '"<<key<<"' in node! The supplied node is not a map (i.e. it is a scalar value, sequence, null, or undefined)";
                     utils_error().raise(LOCAL_INFO,os.str());
                }
                return node[key];
        }

        template <typename... args>
        inline const YAML::Node getVariadicNode(const YAML::Node &node, const std::string &key, const args&... keys)
        {
                if(node.Type()!=YAML::NodeType::Map)
                {
                     std::ostringstream os;
                     os<<"Error accessing YAML key '"<<key<<"' in node! The supplied node is not a map (i.e. it is a scalar value, sequence, null, or undefined)";
                     utils_error().raise(LOCAL_INFO,os.str());
                }
                if(not node[key]) return node[key];
                else return getVariadicNode(node[key], keys...);
        }
        
}

#endif
