//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  Base class for definining the required 
///  contents of Spectrum wrapper objects
///
///  *********************************************
///
///  Authors: 
///  <!-- add name and date if you modify -->
///   
///  \author Ben Farmer
///          (benjamin.farmer@imperial.ac.uk)
///  \date 2016 Feb, 2019 June
///
///  *********************************************

#include <fstream>

#include "gambit/Models/SpectrumContents/spectrum_contents.hpp"
#include "gambit/Models/partmap.hpp"
#include "gambit/Elements/spectrum.hpp"
#include "gambit/Elements/slhaea_helpers.hpp"
#include "gambit/Utils/stream_overloads.hpp"
#include "gambit/Utils/util_functions.hpp"

// Easy name for particle database access
#define PDB Models::ParticleDB()

namespace Gambit 
{
  /// Less confusing name for SLHAea container class
  typedef SLHAea::Coll SLHAstruct;

  namespace SpectrumContents
  {
    /// Function to retrieve all valid indices for this parameter
    /// Returns them in the most general form, vector of vectors, regardless of actual shape
    std::vector<std::vector<int>> Parameter::allowed_indices()
    {
        std::vector<std::vector<int>> out;
        switch(shape().size())
        {
            case 0:
            {
                // No indices, return empty allowed index list
                break; 
            }
            case 1:
            {
                // One index
                for(int i=1; i<=shape().at(0); i++)
                {
                    std::vector<int> indices;
                    indices.push_back(i);
                    out.push_back(indices);
                }
                break;
            }
            case 2:
            {
                // Two indices
                for(int i=1; i<=shape().at(0); i++)
                {
                    for(int j=1; j<=shape().at(1); j++)
                    {
                       std::vector<int> indices;
                       indices.push_back(i);
                       indices.push_back(j);
                       out.push_back(indices);
                    }
                }
                break;
            }
        }
        return out;   
    }
 
    /// Add a parameter to the Contents object
    void Contents::addParameter(const Par::Tags tag, const std::string& name, const std::vector<int>& shape,
                                           const std::string& blockname, const int blockindex)
    {
       parameter_key thiskey(tag,name);
       parameters.emplace(thiskey,Parameter(tag,name,shape,blockname,blockindex));
       /// Special additions for Par::Pole_Mass parameters; uncertainties!
       if(tag==Par::Pole_Mass)
       {
          /// Block is always DMASS, and 'blockindex' is ignored because we use PDG codes to look this up.
          Par::Tags tag_upper(Par::Pole_Mass_1srd_high);
          Par::Tags tag_lower(Par::Pole_Mass_1srd_low);
          parameter_key key_upper(tag_upper,name);
          parameter_key key_lower(tag_lower,name);
          parameters.emplace(key_upper,Parameter(tag_upper,name,shape,"DMASS",blockindex));
          parameters.emplace(key_lower,Parameter(tag_lower,name,shape,"DMASS",blockindex));
       }
    }

    /// Import all parameter definitions from another Contents object
    void Contents::addAllFrom(const Contents& other)
    {
       for(auto p : other.all_parameters())
       {
          parameter_key thiskey(p.tag(),p.name());
          parameters.emplace(thiskey,p);
       }
    }
 
    /// Set the name of this Contents object (i.e. the name of the model to which this spectrum data applies) 
    void Contents::setName(const std::string& name)
    {
       my_name = name;
    }

    /// Function to check if a parameter definition exists in this object, identified by tag and string name
    bool Contents::has_parameter(const Par::Tags tag, const std::string& name) const
    {
       bool found=true;
       if(parameters.find(std::make_pair(tag,name))==parameters.end())
       {
          found=false;
       }
       return found;
    }

    /// Function to check if a parameter definition exists in this object, this time also checking the index shape
    bool Contents::has_parameter(const Par::Tags tag, const std::string& name, const std::vector<int>& indices) const
    {
       bool found=true;
       if(has_parameter(tag,name))
       {
          // Now retrieve parameter definition and check indices
          Parameter par(parameters.at(std::make_pair(tag,name)));
          if(par.shape().size()!=indices.size())
          {
             // Nope, shape is wrong
             found = false;
          }
          else
          {
             // Shape is ok, check ranges (these are all indexed starting at 1)
             for(int i = 1; i<par.shape().size(); i++)
             {
                if(indices.at(i)<1 or indices.at(i)>par.shape().at(i))
                {
                   found = false; // At least one index out of bounds
                }
             }
          }
       }
       else
       {
           found=false;
       }
       return found;
    }

    /// Function to get definition information for one parameter, identified by tag and string name
    Parameter Contents::get_parameter(const Par::Tags tag, const std::string& name) const
    {
       if(not has_parameter(tag,name))
       {
          std::ostringstream errmsg;
          errmsg<<"Could not get_parameter "<<name<<" with tag "<<Par::toString.at(tag)<<" from Contents with name "<<getName()<<std::endl;
          errmsg<<"That parameter has not been defined for this Contents!"<<std::endl;
          utils_error().forced_throw(LOCAL_INFO,errmsg.str());           
       }
       return parameters.at(std::make_pair(tag,name));
    }
 
    /// Function to get block and indices in SLHAea file in which requested index item can be found
    std::pair<std::string,std::vector<int>> Contents::get_SLHA_indices(const Par::Tags tag, const std::string& name, std::vector<int> indices) const
    {
       std::vector<int> SLHA_indices;
       std::string block;
       if(has_parameter(tag, name))
       {
          Parameter p(parameters.at(std::make_pair(tag,name)));
          block = Utils::toUpper(p.blockname());
          if(has_parameter(tag, name, indices))
          {
             // Parameter exists and indices are correct; now convert them to SLHA indices
             switch(p.shape().size())
             {
                case 0:
                {
                   // We have to deal with the MASS block as a special case, because the SLHA indices are PDG codes.
                   if(tag==Par::Pole_Mass & block=="MASS")
                   {
                      std::pair<int, int> pdgpair = PDB.pdg_pair(name);
                      // I don't think we need to worry about the context integer here. Long name should be unique
                      // across all contexts, it is the PDG code that can refer to different things in different model contexts.
                      int pdgcode = pdgpair.first;
                      SLHA_indices.push_back(pdgcode);
                   }
                   if(tag==Par::Pole_Mass_1srd_high or tag==Par::Pole_Mass_1srd_low)
                   {
                      // Pole mass uncertainties also have special retrieval rules, since they are added automatically
                      // Stored with two SLHA indices: one PDG code, and one index indicating 'upper' or 'lower' uncertainty. 
                      std::pair<int, int> pdgpair = PDB.pdg_pair(name);
                      int pdgcode = pdgpair.first;
                      SLHA_indices.push_back(pdgcode);
                      if(tag==Par::Pole_Mass_1srd_high) SLHA_indices.push_back(1);
                      if(tag==Par::Pole_Mass_1srd_low) SLHA_indices.push_back(0);
                   }
                   else
                   {
                      // SLHA index should be directly specified in parameter definition
                      SLHA_indices.push_back(p.blockindex());
                   }
                   break;
                }
                case 1:
                {
                   int index = indices.at(0);
                   // We have to deal with the MASS block as a special case, because the SLHA indices are PDG codes.
                   if(tag==Par::Pole_Mass & block=="MASS")
                   {
                      std::pair<int, int> pdgpair = PDB.pdg_pair(name,index);
                      int pdgcode = pdgpair.first;
                      SLHA_indices.push_back(pdgcode);
                   }
                   if(tag==Par::Pole_Mass_1srd_high or tag==Par::Pole_Mass_1srd_low)
                   {
                      // Pole mass uncertainties also have special retrieval rules, since they are added automatically
                      // Stored with two SLHA indices: one PDG code, and one index indicating 'upper' or 'lower' uncertainty. 
                      std::pair<int, int> pdgpair = PDB.pdg_pair(name,index);
                      int pdgcode = pdgpair.first;
                      SLHA_indices.push_back(pdgcode);
                      if(tag==Par::Pole_Mass_1srd_high) SLHA_indices.push_back(1);
                      if(tag==Par::Pole_Mass_1srd_low) SLHA_indices.push_back(0);
                   }
                   else
                   {
                      // All other SLHA indices are assumed to increment sequentially in tandem with the defined parameter index!
                      // This is a pretty strong assumption/requirement! But I think SLHA conforms to it (aside from the MASS block), 
                      // and we should try to stick to it in our extensions.
                      SLHA_indices.push_back(p.blockindex()+index-1); // Note: parameter indices always start at 1
                   }
                   break;
                }
                case 2:
                {
                   int index1 = indices.at(0);
                   int index2 = indices.at(1);
                   // No weird PDG code stuff for matrix-valued parameters.
                   // blockindex can again be used to offset the SLHA indices from the parameter ones, if needed
                   SLHA_indices.push_back(p.blockindex()+index1-1);
                   SLHA_indices.push_back(p.blockindex()+index2-1);  
                }
             }
          }
          else
          {
             std::ostringstream errmsg;
             errmsg<<"Could not get SLHA indices for parameter "<<name<<" (indices="<<indices<<") with tag "<<Par::toString.at(tag)<<" from Contents with name "<<getName()<<std::endl;
             errmsg<<"A parameter with that name and tag *has* been defined for this Contents, but the requested indices are out of the allowed range (allowed shape has been defined as "<<p.shape()<<")."<<std::endl;
             utils_error().forced_throw(LOCAL_INFO,errmsg.str());           
          }
       }
       else
       {
          std::ostringstream errmsg;
          errmsg<<"Could not get SLHA indices for parameter "<<name<<" (indices="<<indices<<") with tag "<<Par::toString.at(tag)<<" from Contents with name "<<getName()<<std::endl;
          errmsg<<"No parameter with that name and tag has been defined for this Contents!"<<std::endl;
          utils_error().forced_throw(LOCAL_INFO,errmsg.str());           
       }
       return std::make_pair(block,SLHA_indices);
    }


    /// Function to retreive all parameters as a vector
    std::vector<Parameter> Contents::all_parameters() const 
    {
      std::vector<Parameter> all_pars;
      for(auto it=parameters.begin(); it!=parameters.end(); ++it)
      {
         all_pars.push_back(it->second);
      }
      return all_pars; 
    }
    
    /// Function to retreive all parameters matching a certain tag
    std::vector<Parameter> Contents::all_parameters_with_tag(Par::Tags tag) const
    {
      std::vector<Parameter> search_result;
      for(auto& kv : parameters)
      {
        if(kv.second.tag() == tag) search_result.push_back(kv.second);        
      }
      return search_result; 
    }
    
    /// Function to retrieve all parameters matching a certain tag and shape
    std::vector<Parameter> Contents::all_parameters_with_tag_and_shape(Par::Tags tag, std::vector<int>& shape) const
    {
      std::vector<Parameter> search_result;
      for(auto& kv : parameters)
      {
        if(kv.second.tag() == tag and kv.second.shape() == shape) search_result.push_back(kv.second);        
      }
      return search_result; 
    }

    /// Function to retrieve all parameters whose blockname is not SMINPUTS, YUKAWA, CKMBLOCK, or empty.
    std::vector<Parameter> Contents::all_BSM_parameters() const
    {
     std::vector<Parameter> search_result;
     for (auto& kv : parameters)
     {
       if(kv.second.blockname() != "SMINPUTS" || kv.second.blockname() != "YUKAWA" || kv.second.blockname() != "CKMBLOCK" || kv.second.blockname() != "")
       {
         search_result.push_back(kv.second);
       }
     }
     return search_result;
    }


    /// Verify that the supplied Spectrum object conforms to the requirements specified by the Contents class
    void Contents::verify_contents(const Spectrum& spec) const
    {
       const std::vector<Parameter> required_parameters = all_parameters();
       
       for(auto& par : required_parameters)
       {
          const Par::Tags        tag   = par.tag();
          const std::string      name  = par.name();
          const std::vector<int> shape = par.shape();

          // Deal with empty shape case
          if(shape.size()==0)
          {
            // ERROR, please use length 1 vector for scalar case
            std::ostringstream errmsg;           
            errmsg << "Error while verifying contents of Spectrum object against Contents object with name \""<<my_name<<"\" !" << std::endl;
            errmsg << "Encountered a required parameter ("<<Par::toString.at(tag)<<", "<<name<<") with shape.size()==0. This is not allowed; if you want this parameter to be considered a scalar, please set the shape to '1', i.e. std::vector<int> shape = initVector(1). Please fix this parameter in the Contents class." << std::endl;
            utils_error().forced_throw(LOCAL_INFO,errmsg.str());
          }
          // Check scalar case
          else if(shape.size()==1 and shape[0]==1)
          {
            if(not spec.has(tag,name))
            {
              // ERROR, Required parameter not found
              std::ostringstream errmsg;           
              errmsg << "Error while verifying contents of Spectrum object against Contents object with name \""<<my_name<<"\" !" << std::endl;
              errmsg << "Required scalar-valued parameter ("<<Par::toString.at(tag)<<", "<<name<<") is not accessible via subspectrum->get(Par::"<<Par::toString.at(tag)<<", \""<<name<<"\"). Please fix the relevant Spectrum wrapper class so that this parameter can be accessed." << std::endl;
              utils_error().forced_throw(LOCAL_INFO,errmsg.str());
            }
          }
          // Check vector case
          else if(shape.size()==1 and shape[0]>1)
          {
            if(shape[0]<0)
            {
              // ERROR, asked for negative length vector
              std::ostringstream errmsg;           
              errmsg << "Error while verifying contents of Spectrum object against Contents object with name \""<<my_name<<"\" !" << std::endl;
              errmsg << "Encountered a vector-valued required parameter ("<<Par::toString.at(tag)<<", "<<name<<") with negative required length ("<<shape[0]<<")! This is invalid; Please fix this parameter in the Contents class so that the required length is a positive number." << std::endl;
              utils_error().forced_throw(LOCAL_INFO,errmsg.str());
            }
            else
            { 
              for(int i = 1; i<=shape[0]; ++i) {
                if(not spec.has(tag,name,i))
                {
                  // ERROR, Required parameter not found
                  std::ostringstream errmsg;           
                  errmsg << "Error while verifying contents of Spectrum object against Contents object with name \""<<my_name<<"\" !" << std::endl;
                  errmsg << "An entry of the required vector-valued parameter ("<<Par::toString.at(tag)<<", "<<name<<") with required length "<<shape[0]<<" is not accessible via subspectrum->get(Par::"<<Par::toString.at(tag)<<", \""<<name<<"\", "<<i<<"). Please fix the relevant Spectrum wrapper class so that this parameter can be accessed. Keep in mind that you may need to override index_offset() to align the expected indices." << std::endl;
                  utils_error().forced_throw(LOCAL_INFO,errmsg.str());
                }            
              }
            }  
          }
          // Check matrix case
          else if(shape.size()==2)
          {
            if(shape[0]<0 or shape[1]<0)
            {
              // ERROR, asked for negative matrix dimensions
              std::ostringstream errmsg;           
              errmsg << "Error while verifying contents of Spectrum object against Contents object with name \""<<my_name<<"\" !" << std::endl;
              errmsg << "Encountered a matrix-valued required parameter ("<<Par::toString.at(tag)<<", "<<name<<") with at least one negative required dimension (dims = ["<<shape[0]<<", "<<shape[1]<<"])! This is invalid; Please fix the shape settings for this parameter in the Contents so that they are positive numbers." << std::endl;
              utils_error().forced_throw(LOCAL_INFO,errmsg.str());
            }
            else
            {
              for(int i = 1; i<=shape[0]; ++i) {
                for(int j = 1; j<=shape[0]; ++j) {
                   if(not spec.has(tag,name,i,j))
                   {
                     // ERROR, Required parameter not found
                     std::ostringstream errmsg;           
                     errmsg << "Error while verifying contents of Spectrum object against Contents object with name \""<<my_name<<"\" !" << std::endl;
                     errmsg << "An entry of the required matrix-valued parameter ("<<Par::toString.at(tag)<<", "<<name<<") with required dimensions ("<<shape[0]<<", "<<shape[1]<<") is not accessible via subspectrum->get(Par::"<<Par::toString.at(tag)<<", \""<<name<<"\", "<<i<<", "<<j<<"). Please fix the relevant Spectrum wrapper class so that this parameter can be accessed. Keep in mind that you may need to override index_offset() to align the expected indices." << std::endl;
                     utils_error().forced_throw(LOCAL_INFO,errmsg.str());
                   }            
                }
              }  
            }
          }
          // Deal with all other cases
          else
          {
            // ERROR invalid shape
            std::ostringstream errmsg;           
            errmsg << "Error while verifying contents of Spectrum object against Contents object with name \""<<my_name<<"\" !" << std::endl;
            errmsg << "The specified shape for the required parameter ("<<Par::toString.at(tag)<<", "<<name<<") is invalid. The length of the shape vector is only permitted to be 1 or 2 (received shape vector was "<<shape<<"). Please fix this parameter entry in the Contents class."<<std::endl;
            utils_error().forced_throw(LOCAL_INFO,errmsg.str());
          }
       }
       // End constructor
    }

    /// Create template SLHA file to match this Contents: mainly used to help Spectrum object creators to know what exactly is required in the SLHAea objects for a given Contents
    void Contents::create_template_SLHA_file(const std::string& filename) const
    {
       SLHAstruct out;
       for(auto& kv : parameters)
       {
           Par::Tags   tag(kv.first.first);
           std::string name(kv.first.second);
           Parameter par(kv.second);

           // Make sure block exists in output
           SLHAea_check_block(out, par.blockname());

           // Get all allowed indices for this parameter
           std::cout<<"Processing parameter: "<<par.name()<<", "<<Par::toString.at(par.tag())<<", "<<par.shape()<<", "<<par.blockname()<<", "<<par.blockindex()<<std::endl;
           std::vector<std::vector<int>> indices(par.allowed_indices());
           // Get SLHA indices for all entries, and add them to SLHA output
           if(indices.size()==0)
           {
              SLHAea_add(out, par.blockname(), -9999, name + " ("+Par::toString.at(tag)+")");
           }
           else
           {
              for(auto index_list : indices)
              {
                  std::stringstream comment;
                  comment << name;
                  if(index_list.size()>0) comment << index_list;
                  comment << " ("<<Par::toString.at(tag)<<")";

                  std::pair<std::string,std::vector<int>> SLHA_loc = get_SLHA_indices(tag,name,index_list);
                  std::string      block        = SLHA_loc.first;
                  std::vector<int> SLHA_indices = SLHA_loc.second;
                  switch(SLHA_indices.size())
                  {
                      case 0:
                      {
                          // Currently not handled; error
                          std::ostringstream errmsg;
                          errmsg<<"Error generating template SLHA file for Contents "<<getName()<<"! Received null index list from get_SLHA_indices function when looking up indices for parameter "<<name<<" of type "<<Par::toString.at(tag)<<" with spectrum indices "<<index_list<<"! Currently we do not handle the case of SLHA items that have no indices. Please file a bug report if you require this feature, or check the Contents definition if you think this is a mistake.";          
                          utils_error().forced_throw(LOCAL_INFO,errmsg.str());
                          break;
                      }
                      case 1:
                      {
                          int index = SLHA_indices.at(0);
                          SLHAea_add(out, block, index, -9999, comment.str()); 
                          break;
                      }
                      case 2:
                      {
                          int index1 = SLHA_indices.at(0);
                          int index2 = SLHA_indices.at(1);
                          SLHAea_add(out, block, index1, index2, -9999, comment.str());
                      }
                  }
              }
           } 
       }
       // Write to file
       std::ofstream ofs(filename);
       ofs << out;
       ofs.close();
    }
  }
}

#undef PDB
