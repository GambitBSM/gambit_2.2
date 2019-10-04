//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  HDF5 interface reaader class retrieve function
///  overloads.  Add a new overload of the _retrieve
///  function in this file if you want to be able
///  to read a new type for postprocessing.
///
///  *********************************************
///
///  Authors (add name and date if you modify):
///
///  \author Ben Farmer
///          (benjamin.farmer@monash.edu.au)
///  \date 2017 Jan
///
///  \author Pat Scott
///          (p.scott@imperial.ac.uk)
///  \date 2017 March
///
///  *********************************************

#include "gambit/Printers/printers/hdf5reader.hpp"
#include "gambit/Printers/printers/hdf5printer.hpp"

namespace Gambit
{
  namespace Printers
  {

     /// @{ Retrieve functions

     /// Templatable retrieve functions
     #define RETRIEVE(TYPE) _retrieve(TYPE& out, const std::string& l, const uint r, const ulong p) \
        { return  _retrieve_template(out,l,0,r,p); }
     bool HDF5Reader::RETRIEVE(int      )
     bool HDF5Reader::RETRIEVE(uint     )
     bool HDF5Reader::RETRIEVE(long     )
     bool HDF5Reader::RETRIEVE(ulong    )
     bool HDF5Reader::RETRIEVE(longlong )
     bool HDF5Reader::RETRIEVE(ulonglong)
     bool HDF5Reader::RETRIEVE(float    )
     bool HDF5Reader::RETRIEVE(double   )
     #undef RETRIEVE

     // Bools can't quite use the template function directly, since there
     // are some issues with bools and MPI/HDF5 types. Easier to just convert
     // the bool to an int first (this is how they are printed in the first place anyway).
     bool HDF5Reader::_retrieve(bool& out, const std::string& l, const uint rank, const ulong pID)
     {
       uint tmp_out;
       bool tmp_ret;
       tmp_ret = _retrieve_template(tmp_out,l,0,rank,pID);
       out = tmp_out;
       return tmp_ret;
     }

     bool HDF5Reader::_retrieve(ModelParameters& out, const std::string& modelname, const uint rank, const ulong pointID)
     {
        bool is_valid = true;
        /// Work out all the output labels which correspond to the input modelname
        bool found_at_least_one(false);

        //std::cout << "Searching for ModelParameters of model '"<<modelname<<"'"<<std::endl;
        // Iterate through names in HDF5 group
        for(std::vector<std::string>::const_iterator
            it = all_datasets.begin();
            it!= all_datasets.end(); ++it)
        {
          //std::cout << "Candidate: " <<*it<<std::endl;
          std::string param_name; // *output* of parsing function, parameter name
          std::string label_root; // *output* of parsing function, label minus parameter name
          if(parse_label_for_ModelParameters(*it, modelname, param_name, label_root))
          {
            // Add the found parameter name to the ModelParameters object
            out._definePar(param_name);
            if(found_at_least_one)
            {
              if(out.getOutputName()!=label_root)
              {
                std::ostringstream err;
                err << "Error! HDF5Reader could not retrieve ModelParameters matching the model name '"
                    <<modelname<<"' in the HDF5 file:group "<<file<<":"<<group
                    <<"' (while calling 'retrieve'). Candidate parameters WERE found, however their dataset "
                    <<"labels indicate the presence of an inconsistency or ambiguity in the output. For "
                    <<"example, we just tried to retrive a model parameter from the dataset:\n  "<<*it
                    <<"\nand successfully found the parameter "<<param_name
                    <<", however the root of the label, that is,\n  "<<label_root
                    <<"\ndoes not match the root expected based upon previous parameter retrievals for this "
                    <<"model, which was\n  "<<out.getOutputName()<<"\nThis may indicate that multiple sets "
                    <<"of model parameters are present in the output file for the same model! This is not "
                    <<"allowed, please report this bug against whatever master YAML file (or external code?) "
                    <<"produced the output file you are trying to read.";
                printer_error().raise(LOCAL_INFO,err.str());
              }
            }
            else
            {
              out.setOutputName(label_root);
            }
            // Get the corresponding value out of the data file
            double value; // *output* of retrieve function
            bool tmp_is_valid;
            tmp_is_valid = _retrieve(value, *it, rank, pointID);
            found_at_least_one = true;
            if(tmp_is_valid)
            {
               out.setValue(param_name, value);
            }
            else
            {
               // If one parameter value is 'invalid' then we cannot reconstruct
               // the ModelParameters object, so we mark the whole thing invalid.
               out.setValue(param_name, 0);
               is_valid = false;
            }
          }
        }

        if(not found_at_least_one)
        {
          // Didn't find any matches!
           std::ostringstream err;
           err << "Error! HDF5Reader failed to find any ModelParameters matching the model name '"<<modelname<<"' in the HDF5 file:group "<<file<<":"<<group<<"' (while calling 'retrieve'). Please check that model name and input file/group are correct.";
           printer_error().raise(LOCAL_INFO,err.str());
        }
        /// done!
        return is_valid;
     }

     /// Retrieve MSSM spectrum information as an SLHAea object
     bool HDF5Reader::_retrieve(MSSM_SLHAstruct& out, const std::string& spec_type, const uint rank, const ulong pointID)
     {
        bool is_valid = true;
        /// Work out all the output labels which correspond to the input modelname
        bool found_at_least_one(false);

        //std::cout << "Searching for Spectrum entries for spectrum '"<<spec_type<<"'"<<std::endl;
        // Iterate through names in HDF5 group
        for(std::vector<std::string>::const_iterator
            it = all_datasets.begin();
            it!= all_datasets.end(); ++it)
        {
          //std::cout << "Candidate: " <<*it<<std::endl;
          std::string entry; // *output* of parsing function, spectrum entry name
          std::string tag;   // *output* of parsing function, spectrum entry tag
          std::string label_root; // *output* of parsing function, label minus entry name and tag
                  
          if(parse_label_for_spectrum_entry(*it, spec_type, "SpecBit", "get_MSSM_spectrum_as_map", entry, tag, label_root))
          {
            // Ok, found! Now retrieve the data
            double value; // *output* of retrieve function
            bool tmp_is_valid;
            tmp_is_valid = _retrieve(value, *it, rank, pointID);
            found_at_least_one = true;
            if(tmp_is_valid)
            {
              // Figure out where in the SLHAea object we need to put this entry.
              // TODO: After the SpecBit redesign this should be possible to automate, since there will be routines that specify
              // exactly what blocks certain spectrum information should live in. But for now we'll have to do it
              // manually with a massive 'if' block
              if(Utils::iequals(tag,"Pole_Mass"))
              {
                 SLHAea_check_block(out, "MASS");
                 if     (Utils::iequals(entry, "A0")) SLHAea_add(out, "MASS", 36, value, entry+" ("+tag+")");
                 else if(Utils::iequals(entry, "H+")) SLHAea_add(out, "MASS", 37, value, entry+" ("+tag+")");
                 else if(Utils::iequals(entry, "W+")) SLHAea_add(out, "MASS", 24, value, entry+" ("+tag+")");
                 else if(Utils::iequals(entry, "~g")) SLHAea_add(out, "MASS", 1000021, value, entry+" ("+tag+")");
                 else if(Utils::iequals(entry, "Z0")) ; // Do nothing, get this from SMINPUTS
                 else if(Utils::iequals(entry, "b")) ;
                 else if(Utils::iequals(entry, "d_3")) ;
                 else if(Utils::iequals(entry, "u_3")) ;
                 else if(Utils::iequals(entry, "e-")) ;
                 else if(Utils::iequals(entry, "mu-")) ;
                 else if(Utils::iequals(entry, "tau-")) ;
                 else if(Utils::iequals(entry, "e-_1")) ;
                 else if(Utils::iequals(entry, "e-_2")) ;
                 else if(Utils::iequals(entry, "e-_3")) ;
                 else if(Utils::iequals(entry, "nu_1")) ;
                 else if(Utils::iequals(entry, "nu_2")) ;
                 else if(Utils::iequals(entry, "nu_3")) ;
                 else if(Utils::iequals(entry, "g")) ;
                 else if(Utils::iequals(entry, "gamma")) ;
                 else if(Utils::iequals(entry, "t")) ;
                 else if(Utils::iequals(entry, "g")) ;
                 else if(Utils::iequals(entry, "g")) ;
                 else if(Utils::iequals(entry, "h0_1")) SLHAea_add(out, "MASS", 25, value, entry+" ("+tag+")");
                 else if(Utils::iequals(entry, "h0_2")) SLHAea_add(out, "MASS", 35, value, entry+" ("+tag+")");
                 else if(Utils::iequals(entry, "~chi+_1")) SLHAea_add(out, "MASS", 1000024, value, entry+" ("+tag+")");
                 else if(Utils::iequals(entry, "~chi+_2")) SLHAea_add(out, "MASS", 1000037, value, entry+" ("+tag+")");
                 else if(Utils::iequals(entry, "~chi0_1")) SLHAea_add(out, "MASS", 1000022, value, entry+" ("+tag+")");
                 else if(Utils::iequals(entry, "~chi0_2")) SLHAea_add(out, "MASS", 1000023, value, entry+" ("+tag+")");
                 else if(Utils::iequals(entry, "~chi0_3")) SLHAea_add(out, "MASS", 1000025, value, entry+" ("+tag+")");
                 else if(Utils::iequals(entry, "~chi0_4")) SLHAea_add(out, "MASS", 1000035, value, entry+" ("+tag+")");
                 else if(Utils::iequals(entry, "~d_1")) SLHAea_add(out, "MASS", 1000001, value, entry+" ("+tag+")");
                 else if(Utils::iequals(entry, "~d_2")) SLHAea_add(out, "MASS", 1000003, value, entry+" ("+tag+")");
                 else if(Utils::iequals(entry, "~d_3")) SLHAea_add(out, "MASS", 1000005, value, entry+" ("+tag+")");
                 else if(Utils::iequals(entry, "~d_4")) SLHAea_add(out, "MASS", 2000001, value, entry+" ("+tag+")");
                 else if(Utils::iequals(entry, "~d_5")) SLHAea_add(out, "MASS", 2000003, value, entry+" ("+tag+")");
                 else if(Utils::iequals(entry, "~d_6")) SLHAea_add(out, "MASS", 2000005, value, entry+" ("+tag+")");
                 else if(Utils::iequals(entry, "~u_1")) SLHAea_add(out, "MASS", 1000002, value, entry+" ("+tag+")");
                 else if(Utils::iequals(entry, "~u_2")) SLHAea_add(out, "MASS", 1000004, value, entry+" ("+tag+")");
                 else if(Utils::iequals(entry, "~u_3")) SLHAea_add(out, "MASS", 1000006, value, entry+" ("+tag+")");
                 else if(Utils::iequals(entry, "~u_4")) SLHAea_add(out, "MASS", 2000002, value, entry+" ("+tag+")");
                 else if(Utils::iequals(entry, "~u_5")) SLHAea_add(out, "MASS", 2000004, value, entry+" ("+tag+")");
                 else if(Utils::iequals(entry, "~u_6")) SLHAea_add(out, "MASS", 2000006, value, entry+" ("+tag+")");
                 else if(Utils::iequals(entry, "~e-_1")) SLHAea_add(out, "MASS", 1000011, value, entry+" ("+tag+")");
                 else if(Utils::iequals(entry, "~e-_2")) SLHAea_add(out, "MASS", 1000013, value, entry+" ("+tag+")");
                 else if(Utils::iequals(entry, "~e-_3")) SLHAea_add(out, "MASS", 1000015, value, entry+" ("+tag+")");
                 else if(Utils::iequals(entry, "~e-_4")) SLHAea_add(out, "MASS", 2000011, value, entry+" ("+tag+")");
                 else if(Utils::iequals(entry, "~e-_5")) SLHAea_add(out, "MASS", 2000013, value, entry+" ("+tag+")");
                 else if(Utils::iequals(entry, "~e-_6")) SLHAea_add(out, "MASS", 2000015, value, entry+" ("+tag+")");
                 else if(Utils::iequals(entry, "~nu_1")) SLHAea_add(out, "MASS", 1000012, value, entry+" ("+tag+")");
                 else if(Utils::iequals(entry, "~nu_2")) SLHAea_add(out, "MASS", 1000014, value, entry+" ("+tag+")");
                 else if(Utils::iequals(entry, "~nu_3")) SLHAea_add(out, "MASS", 1000016, value, entry+" ("+tag+")");
                 else
                 {
                    std::ostringstream err;
                    err << "Error! HDF5Reader encountered an error while attempting to read a spectrum of type '"<<spec_type<<"' from the HDF5 file:group "<<file<<":"<<group<<"' (while calling 'retrieve'). An unrecognised entry with tag 'Pole_Mass' was encountered (label was "<<entry<<")."<<std::endl<<"The full HDF5 entry label was: "<<*it;
                    printer_error().raise(LOCAL_INFO,err.str());
                 }
              }
              else if(Utils::iequals(tag,"mass1"))
              {
                 SLHAea_check_block(out, "MSOFT");
                 SLHAea_check_block(out, "HMIX");
                 if     (Utils::iequals(entry, "M1")) SLHAea_add(out, "MSOFT", 1, value, entry+" ("+tag+")");
                 else if(Utils::iequals(entry, "M2")) SLHAea_add(out, "MSOFT", 2, value, entry+" ("+tag+")");
                 else if(Utils::iequals(entry, "M3")) SLHAea_add(out, "MSOFT", 3, value, entry+" ("+tag+")");
                 else if(Utils::iequals(entry, "Mu")) SLHAea_add(out, "HMIX", 1, value, entry+" ("+tag+")");
                 else if(Utils::iequals(entry, "d_1")) ;
                 else if(Utils::iequals(entry, "d_2")) ;
                 else if(Utils::iequals(entry, "d_3")) ;
                 else if(Utils::iequals(entry, "u_1")) ;
                 else if(Utils::iequals(entry, "u_2")) ;
                 else if(Utils::iequals(entry, "u_3")) ;
                 else if(Utils::iequals(entry, "e-_1")) ;
                 else if(Utils::iequals(entry, "e-_2")) ;
                 else if(Utils::iequals(entry, "e-_3")) ;
                 else if(Utils::iequals(entry, "g")) ;
                 else if(Utils::iequals(entry, "gamma")) ;
                 //...etc. TODO: incomplete!
              } 
            }
            else
            {
              // If one entry is 'invalid' then we cannot reconstruct
              // the Spectrum object, so we mark the whole thing invalid.
              is_valid = false;
            }
          }
        }

        if(not found_at_least_one)
        {
          // Didn't find any matches!
           std::ostringstream err;
           err << "Error! HDF5Reader failed to find any Spectrum entries matching the spectrum type '"<<spec_type<<"' in the HDF5 file:group "<<file<<":"<<group<<"' (while calling 'retrieve'). Please check that spectrum type and input file/group are correct.";
           printer_error().raise(LOCAL_INFO,err.str());
        }
        /// done!
        return is_valid;
     }


     bool HDF5Reader::_retrieve(std::vector<double>& /*out*/,  const std::string& /*label*/, const uint /*rank*/, const ulong /*pointID*/)
     { printer_error().raise(LOCAL_INFO,"NOT YET IMPLEMENTED"); return false; }
     bool HDF5Reader::_retrieve(map_str_dbl& /*out*/,          const std::string& /*label*/, const uint /*rank*/, const ulong /*pointID*/)
     { printer_error().raise(LOCAL_INFO,"NOT YET IMPLEMENTED"); return false; }
     bool HDF5Reader::_retrieve(triplet<double>& /*out*/,      const std::string& /*label*/, const uint /*rank*/, const ulong /*pointID*/)
     { printer_error().raise(LOCAL_INFO,"NOT YET IMPLEMENTED"); return false; }
     bool HDF5Reader::_retrieve(map_intpair_dbl& /*out*/,      const std::string& /*label*/, const uint /*rank*/, const ulong /*pointID*/)
     { printer_error().raise(LOCAL_INFO,"NOT YET IMPLEMENTED"); return false; }

    #ifndef SCANNER_STANDALONE // All the types inside HDF5_MODULE_BACKEND_TYPES need to go inside this def guard.

       bool HDF5Reader::_retrieve(DM_nucleon_couplings& /*out*/, const std::string& /*label*/, const uint /*rank*/, const ulong /*pointID*/)
       { printer_error().raise(LOCAL_INFO,"NOT YET IMPLEMENTED"); return false; }
       bool HDF5Reader::_retrieve(Flav_KstarMuMu_obs& /*out*/, const std::string& /*label*/, const uint /*rank*/, const ulong /*pointID*/)
       { printer_error().raise(LOCAL_INFO,"NOT YET IMPLEMENTED"); return false; }

     #endif

     /// @}

  }
}
