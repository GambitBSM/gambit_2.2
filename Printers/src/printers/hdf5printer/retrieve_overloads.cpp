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

     struct SLHAcombo
     {
        SLHAcombo(const std::string& t, const std::string& b, int i)
         : tag(t)
         , block(b)
         , indices{i}
        {}

        SLHAcombo(const std::string& t, const std::string& b, int i, int j)
         : tag(t)
         , block(b)
         , indices{i,j}
        {}

        SLHAcombo() : tag(), block(), indices() {}

        std::string tag;
        std::string block;
        std::vector<int> indices;
     };

     bool HDF5Reader::retrieve_and_add_to_SLHAea(SLHAstruct& out, bool& found, const std::string& spec_type, const std::string& entry, const SLHAcombo& item, const std::set<std::string>& all_dataset_labels, const uint rank, const ulong pointID)
     {
        std::string tag   = item.tag;
        std::string block = item.block;
        std::vector<int> indices = item.indices;

        // Create full dataset label
        std::stringstream dataset_label;
        dataset_label<<"#"<<spec_type<<" @SpecBit::get_MSSM_spectrum_as_map::"<<entry;
        if(tag!="") dataset_label<<" "<<tag;

        std::cout<<dataset_label.str()<<std::endl;
        auto jt = all_dataset_labels.find(dataset_label.str());
        if(jt==all_dataset_labels.end())
        {
           found = false; // No entry with this name!
           return false;
        }
        else
        {
           found = true;
        }

        // Ok, found! Now retrieve the data
        double value = -999; // *output* of retrieve function
        bool tmp_is_valid = false;
        tmp_is_valid = _retrieve(value, dataset_label.str(), rank, pointID);
        std::cout<<"Spectrum entry found! entry:"<<entry<<", tag:"<<tag<<", valid:"<<tmp_is_valid<<", value:"<<value<<std::endl;
        if(tmp_is_valid)
        {
            // Stick entry into the SLHAea object
            SLHAea_check_block(out, block); // Make sure block exists first
            if(indices.size()==1)
            { 
                SLHAea_add(out, block, indices.at(0), value, entry+" ("+tag+")");
            }
            else if(indices.size()==2)
            {
                SLHAea_add(out, block, indices.at(0), indices.at(1), value, entry+" ("+tag+")");
            }
            else
            {
                std::ostringstream err;
                err<<"Received invalid number of target SLHA indices for dataset: "<<dataset_label.str()<<std::endl<<"Indices were: "<<indices;
                printer_error().raise(LOCAL_INFO,err.str()); 
            }
        }
        return tmp_is_valid; 
     }

     /// Retrieve MSSM spectrum information as an SLHAea object
     bool HDF5Reader::_retrieve(MSSM_SLHAstruct& out_main, const std::string& spec_type, const uint rank, const ulong pointID)
     {
        SLHAstruct& out(out_main); // Interpret as ordinary SLHAea base class to get operator[] etc
        bool is_valid = true;

        // Rather than iterate through the datasets, we know what entries we need to find, so we will just
        // directly look for them. 
        // TODO: We can automate this after the SpecBit redesign, and probably
        // just use the spectrum "setter" functions to insert this data directly into Spectrum objects.
        // Unfortunately those don't exist in the current SimpleSpectrum objects, but they will exist after
        // the redesign.
        std::map<std::string,SLHAcombo> labels_to_SLHA;
 
        // MASS
        labels_to_SLHA["A0"     ] = SLHAcombo("Pole_Mass", "MASS", 36);
        labels_to_SLHA["H+"     ] = SLHAcombo("Pole_Mass", "MASS", 37);
        labels_to_SLHA["W+"     ] = SLHAcombo("Pole_Mass", "MASS", 24);
        labels_to_SLHA["h0_1"   ] = SLHAcombo("Pole_Mass", "MASS", 25);
        labels_to_SLHA["h0_2"   ] = SLHAcombo("Pole_Mass", "MASS", 35);
        labels_to_SLHA["~g"     ] = SLHAcombo("Pole_Mass", "MASS", 1000021);
        labels_to_SLHA["~chi+_1"] = SLHAcombo("Pole_Mass", "MASS", 1000024);
        labels_to_SLHA["~chi+_2"] = SLHAcombo("Pole_Mass", "MASS", 1000037);
        labels_to_SLHA["~chi0_1"] = SLHAcombo("Pole_Mass", "MASS", 1000022);
        labels_to_SLHA["~chi0_2"] = SLHAcombo("Pole_Mass", "MASS", 1000023);
        labels_to_SLHA["~chi0_3"] = SLHAcombo("Pole_Mass", "MASS", 1000025);
        labels_to_SLHA["~chi0_4"] = SLHAcombo("Pole_Mass", "MASS", 1000035);
        labels_to_SLHA["~d_1"   ] = SLHAcombo("Pole_Mass", "MASS", 1000001);
        labels_to_SLHA["~d_2"   ] = SLHAcombo("Pole_Mass", "MASS", 1000003);
        labels_to_SLHA["~d_3"   ] = SLHAcombo("Pole_Mass", "MASS", 1000005);
        labels_to_SLHA["~d_4"   ] = SLHAcombo("Pole_Mass", "MASS", 2000001);
        labels_to_SLHA["~d_5"   ] = SLHAcombo("Pole_Mass", "MASS", 2000003);
        labels_to_SLHA["~d_6"   ] = SLHAcombo("Pole_Mass", "MASS", 2000005);
        labels_to_SLHA["~u_1"   ] = SLHAcombo("Pole_Mass", "MASS", 1000002);
        labels_to_SLHA["~u_2"   ] = SLHAcombo("Pole_Mass", "MASS", 1000004);
        labels_to_SLHA["~u_3"   ] = SLHAcombo("Pole_Mass", "MASS", 1000006);
        labels_to_SLHA["~u_4"   ] = SLHAcombo("Pole_Mass", "MASS", 2000002);
        labels_to_SLHA["~u_5"   ] = SLHAcombo("Pole_Mass", "MASS", 2000004);
        labels_to_SLHA["~u_6"   ] = SLHAcombo("Pole_Mass", "MASS", 2000006);
        labels_to_SLHA["~e-_1"  ] = SLHAcombo("Pole_Mass", "MASS", 1000011);
        labels_to_SLHA["~e-_2"  ] = SLHAcombo("Pole_Mass", "MASS", 1000013);
        labels_to_SLHA["~e-_3"  ] = SLHAcombo("Pole_Mass", "MASS", 1000015);
        labels_to_SLHA["~e-_4"  ] = SLHAcombo("Pole_Mass", "MASS", 2000011);
        labels_to_SLHA["~e-_5"  ] = SLHAcombo("Pole_Mass", "MASS", 2000013);
        labels_to_SLHA["~e-_6"  ] = SLHAcombo("Pole_Mass", "MASS", 2000015);
        labels_to_SLHA["~nu_1"  ] = SLHAcombo("Pole_Mass", "MASS", 1000012);
        labels_to_SLHA["~nu_2"  ] = SLHAcombo("Pole_Mass", "MASS", 1000014);
        labels_to_SLHA["~nu_3"  ] = SLHAcombo("Pole_Mass", "MASS", 1000016);

        // MSOFT
        labels_to_SLHA["M1"  ] = SLHAcombo("mass1", "MSOFT", 1);
        labels_to_SLHA["M2"  ] = SLHAcombo("mass1", "MSOFT", 2);
        labels_to_SLHA["M3"  ] = SLHAcombo("mass1", "MSOFT", 3);
        labels_to_SLHA["mHd2"] = SLHAcombo("mass2", "MSOFT", 21);
        labels_to_SLHA["mHu2"] = SLHAcombo("mass2", "MSOFT", 22);
    
        // HMIX
        labels_to_SLHA["Mu"]  = SLHAcombo("mass1", "HMIX", 1);
        // Need these two for Higgs vev and tanbeta. Not SLHA, so store in TEMP block temporarily.
        //labels_to_SLHA["vd"]  = SLHAcombo("mass1", "TEMP", 1);
        //labels_to_SLHA["vu"]  = SLHAcombo("mass1", "TEMP", 2);
        //labels_to_SLHA["mA2"] = SLHAcombo("mass2", "HMIX", 4);
 
        // TD, TU, TE 
        #define LABELNXN(N,baseentry,tag,block) \
          for(int i=1; i<=N; i++){ for(int j=1; j<=N; j++) { \
            std::stringstream entry; \
            entry<<baseentry<<"_("<<i<<","<<j<<")"; \
            labels_to_SLHA[entry.str()] = SLHAcombo(tag, block, i, j); \
          }}
        LABELNXN(3,"TYd","mass1","TD")
        LABELNXN(3,"TYu","mass1","TU")
        LABELNXN(3,"TYe","mass1","TE")

        // MSQ2, MSL2, MSD2, MSU2, MSE2
        LABELNXN(3,"mq2","mass2","MSQ2")
        LABELNXN(3,"ml2","mass2","MSL2")
        LABELNXN(3,"md2","mass2","MSD2")
        LABELNXN(3,"mu2","mass2","MSU2")
        LABELNXN(3,"me2","mass2","MSE2")

        // NMIX, UMIX, VMIX
        LABELNXN(4,"~chi0","Pole_Mixing","NMIX")
        LABELNXN(2,"~chi-","Pole_Mixing","UMIX")
        LABELNXN(2,"~chi+","Pole_Mixing","VMIX")
    
        // USQMIX, DSQMIX, SELMIX, SNUMIX
        LABELNXN(6,"~u","Pole_Mixing","USQMIX")
        LABELNXN(6,"~d","Pole_Mixing","DSQMIX")
        LABELNXN(6,"~e-","Pole_Mixing","SELMIX")
        LABELNXN(3,"~nu","Pole_Mixing","SNUMIX")
        #undef LABELNXN

        // YD, YU, YE
        #define LABEL3X3DIAG(baseentry,tag,block) \
          for(int i=1; i<=3; i++){ \
            std::stringstream entry; \
            entry<<baseentry<<"_("<<i<<","<<i<<")"; \
            labels_to_SLHA[entry.str()] = SLHAcombo(tag, block, i, i); \
          }
        LABEL3X3DIAG("Yd","dimensionless","YD")
        LABEL3X3DIAG("Yu","dimensionless","YU")
        LABEL3X3DIAG("Ye","dimensionless","YE")
        #undef LABEL3X3DIAG

        // GAUGE
        labels_to_SLHA["g1"] = SLHAcombo("dimensionless", "TEMP", 3); // Convert to g'
        labels_to_SLHA["g2"] = SLHAcombo("dimensionless", "GAUGE", 2);
        labels_to_SLHA["g3"] = SLHAcombo("dimensionless", "GAUGE", 3);

        // Read all dataset labels into a structure that we can search quickly
        std::set<std::string> all_dataset_labels = get_all_labels();

        // Read the "scale" entry first, since we need to add this info to the block
        // top rows.
 
        bool found(true);
        retrieve_and_add_to_SLHAea(out, found, spec_type, "scale(Q)", SLHAcombo("", "TEMP", 0), all_dataset_labels, rank, pointID);
        double scale;
        if(not found)
        {
           // In some older datasets we forgot to add the scale to the output.
           // For now we will assume the spectrum was output by FlexibleSUSY, in which case
           // the running parameters will be defined at the SUSY scale (geometric mean of
           // DRbar stop masses). TODO: Set this behvaiour with an option, maybe? Not sure how though.
           
           // TODO: This should come from the DRbar squark mass matrices, but I need to check with Peter
           // how it should be calculated. For now I will used the stop pole masses, for testing. TODO: actually
           // I cannot even easily tell which ones are the stop pole masses. They are probably the lightest two though,
           // so I will just assume that for now.

           double ms1,ms2=1e99;
           for(int i=1; i<=6; i++)
           {
               std::stringstream ss;
               ss<<"~u_"<<i;
               bool found_temp(true);
               retrieve_and_add_to_SLHAea(out, found_temp, spec_type, ss.str(), SLHAcombo("Pole_Mass", "TEMP", 10+i), all_dataset_labels, rank, pointID);
               if(not found_temp)
               {
                   std::ostringstream err;
                   err<<"Failed to find "<<ss.str()<<" (Pole_Mass) need for scale calculation!";
                   printer_error().raise(LOCAL_INFO,err.str());
               }
               double m = SLHAea_get(out,"TEMP",10+i);
               if(m<ms1)
               {
                   ms2 = ms1;
                   ms1 = m;
               } 
               else if(m<ms2)
               {
                   ms2 = m;
               }
           }
           scale = sqrt(ms1*ms2);   
        } 
        else
        { 
           scale = SLHAea_get(out,"TEMP",0);
        } 

        // Add blocks that require scale info
        SLHAea_add_block(out, "GAUGE", scale);
        SLHAea_add_block(out, "YU", scale);
        SLHAea_add_block(out, "YD", scale);
        SLHAea_add_block(out, "YE", scale);
        SLHAea_add_block(out, "TU", scale);
        SLHAea_add_block(out, "TD", scale);
        SLHAea_add_block(out, "TE", scale);
        SLHAea_add_block(out, "HMIX", scale);
        SLHAea_add_block(out, "MSQ2", scale);
        SLHAea_add_block(out, "MSL2", scale);
        SLHAea_add_block(out, "MSD2", scale);
        SLHAea_add_block(out, "MSU2", scale);
        SLHAea_add_block(out, "MSE2", scale);
        SLHAea_add_block(out, "MSOFT", scale);

        // Another manual intervention is required for mA2. Seems we didn't used to save this.
        // But again, we can compute it from other stuff: m^2_A=[m3^2/cosBsinB],
        // where m3^2 is what we call BMu (or small b).
        // Need tanb for this too.
        #define GETPAR(OUT,NAME,TAG,TEMP_INDEX) \
        { \
           bool found_tmp; \
           retrieve_and_add_to_SLHAea(out, found_tmp, spec_type, NAME, SLHAcombo(TAG, "TEMP", TEMP_INDEX), all_dataset_labels, rank, pointID); \
           if(not found_tmp) \
           { \
              std::ostringstream err; \
              err<<"Failed to find "<<NAME<<" ("<<TAG<<") need to compute HMIX block entries!"; \
              printer_error().raise(LOCAL_INFO,err.str()); \
           } \
           OUT = SLHAea_get(out,"TEMP",TEMP_INDEX); \
        }

        double vd,vu,BMu;
        GETPAR(vd,"vd","mass1",1)
        GETPAR(vu,"vu","mass1",2)
        GETPAR(BMu,"BMu","mass2",4)
        #undef GETPAR

        double tb = vu/vd;
        double cb = cos(atan(tb));
        double sb = sin(atan(tb));
        double vev = sqrt(vu*vu + vd*vd);
        double mA2 = BMu / (cb*sb);
        
        // Automatically extract and add the rest of the entries
        for(auto it=labels_to_SLHA.begin(); it!=labels_to_SLHA.end(); ++it)
        {
           bool found(true);
           bool tmp_is_valid = retrieve_and_add_to_SLHAea(out, found, spec_type, it->first, it->second, all_dataset_labels, rank, pointID);
           if(not found)
           {
              std::ostringstream err;
              err << "Error! HDF5Reader encountered an error while attempting to read a spectrum of type '"<<spec_type<<"' from the HDF5 file:group "<<file<<":"<<group<<"' (while calling 'retrieve'). A required dataset could not be found ("<<it->first<<")";
              printer_error().raise(LOCAL_INFO,err.str());
           }
           else if(not tmp_is_valid)
           {
              // No need to read any more if some required spectrum entries are invalid. Whole spectrum is invalid.
              is_valid = false;
              break;
           }
        }

        // Need to manually fix up a few entries where we didn't store the spectrum info directly
        // in SLHA format.
        if(is_valid)
        {
           SLHAea_add(out, "HMIX", 2, tb, "tan beta (Q)");
           SLHAea_add(out, "HMIX", 3, vev, "Higgs vev (Q)");
           SLHAea_add(out, "HMIX", 4, mA2, "m_A^2 = BMu/(cb*sb) (Q)");
           // Normalisation of g1
           double gprime = SLHAea_get(out,"TEMP",3)*sqrt(3./5.);
           SLHAea_add(out, "GAUGE", 1, gprime, "g' (Q)", true);
        }

        // //std::cout << "Searching for Spectrum entries for spectrum '"<<spec_type<<"'"<<std::endl;
        // // Iterate through names in HDF5 group
        // for(std::vector<std::string>::const_iterator
        //     it = all_datasets.begin();
        //     it!= all_datasets.end(); ++it)
        // {
        //   std::cout << "Candidate: " <<*it<<std::endl;
        //   std::string entry; // *output* of parsing function, spectrum entry name
        //   std::string tag;   // *output* of parsing function, spectrum entry tag
        //   std::string label_root; // *output* of parsing function, label minus entry name and tag
        //           
        //   if(parse_label_for_spectrum_entry(*it, spec_type, "SpecBit", "get_MSSM_spectrum_as_map", entry, tag, label_root))
        //   {
        //     // Ok, found! Now retrieve the data
        //     double value = -999; // *output* of retrieve function
        //     bool tmp_is_valid = false;
        //     tmp_is_valid = _retrieve(value, *it, rank, pointID);
        //     found_at_least_one = true;
        //     std::cout<<"Spectrum entry found! entry:"<<entry<<", tag:"<<tag<<", valid:"<<tmp_is_valid<<", value:"<<value<<std::endl;
        //     if(tmp_is_valid)
        //     {
        //       // Figure out where in the SLHAea object we need to put this entry.
        //       // TODO: After the SpecBit redesign this should be possible to automate, since there will be routines that specify
        //       // exactly what blocks certain spectrum information should live in. But for now we'll have to do it
        //       // manually with a massive 'if' block
        //       if(Utils::iequals(tag,"Pole_Mass"))
        //       {
        //          SLHAea_check_block(out, "MASS");
        //          if     (Utils::iequals(entry, "A0")) SLHAea_add(out, "MASS", 36, value, entry+" ("+tag+")");
        //          else if(Utils::iequals(entry, "H+")) SLHAea_add(out, "MASS", 37, value, entry+" ("+tag+")");
        //          else if(Utils::iequals(entry, "W+")) SLHAea_add(out, "MASS", 24, value, entry+" ("+tag+")");
        //          else if(Utils::iequals(entry, "~g")) SLHAea_add(out, "MASS", 1000021, value, entry+" ("+tag+")");
        //          else if(Utils::iequals(entry, "Z0")) ; // Do nothing, get this from SMINPUTS
        //          else if(Utils::iequals(entry, "b")) ;
        //          else if(Utils::iequals(entry, "d_3")) ;
        //          else if(Utils::iequals(entry, "u_3")) ;
        //          else if(Utils::iequals(entry, "e-")) ;
        //          else if(Utils::iequals(entry, "mu-")) ;
        //          else if(Utils::iequals(entry, "tau-")) ;
        //          else if(Utils::iequals(entry, "e-_1")) ;
        //          else if(Utils::iequals(entry, "e-_2")) ;
        //          else if(Utils::iequals(entry, "e-_3")) ;
        //          else if(Utils::iequals(entry, "nu_1")) ;
        //          else if(Utils::iequals(entry, "nu_2")) ;
        //          else if(Utils::iequals(entry, "nu_3")) ;
        //          else if(Utils::iequals(entry, "g")) ;
        //          else if(Utils::iequals(entry, "gamma")) ;
        //          else if(Utils::iequals(entry, "t")) ;
        //          else if(Utils::iequals(entry, "g")) ;
        //          else if(Utils::iequals(entry, "g")) ;
        //          else if(Utils::iequals(entry, "h0_1")) SLHAea_add(out, "MASS", 25, value, entry+" ("+tag+")");
        //          else if(Utils::iequals(entry, "h0_2")) SLHAea_add(out, "MASS", 35, value, entry+" ("+tag+")");
        //          else if(Utils::iequals(entry, "~chi+_1")) SLHAea_add(out, "MASS", 1000024, value, entry+" ("+tag+")");
        //          else if(Utils::iequals(entry, "~chi+_2")) SLHAea_add(out, "MASS", 1000037, value, entry+" ("+tag+")");
        //          else if(Utils::iequals(entry, "~chi0_1")) SLHAea_add(out, "MASS", 1000022, value, entry+" ("+tag+")");
        //          else if(Utils::iequals(entry, "~chi0_2")) SLHAea_add(out, "MASS", 1000023, value, entry+" ("+tag+")");
        //          else if(Utils::iequals(entry, "~chi0_3")) SLHAea_add(out, "MASS", 1000025, value, entry+" ("+tag+")");
        //          else if(Utils::iequals(entry, "~chi0_4")) SLHAea_add(out, "MASS", 1000035, value, entry+" ("+tag+")");
        //          else if(Utils::iequals(entry, "~d_1")) SLHAea_add(out, "MASS", 1000001, value, entry+" ("+tag+")");
        //          else if(Utils::iequals(entry, "~d_2")) SLHAea_add(out, "MASS", 1000003, value, entry+" ("+tag+")");
        //          else if(Utils::iequals(entry, "~d_3")) SLHAea_add(out, "MASS", 1000005, value, entry+" ("+tag+")");
        //          else if(Utils::iequals(entry, "~d_4")) SLHAea_add(out, "MASS", 2000001, value, entry+" ("+tag+")");
        //          else if(Utils::iequals(entry, "~d_5")) SLHAea_add(out, "MASS", 2000003, value, entry+" ("+tag+")");
        //          else if(Utils::iequals(entry, "~d_6")) SLHAea_add(out, "MASS", 2000005, value, entry+" ("+tag+")");
        //          else if(Utils::iequals(entry, "~u_1")) SLHAea_add(out, "MASS", 1000002, value, entry+" ("+tag+")");
        //          else if(Utils::iequals(entry, "~u_2")) SLHAea_add(out, "MASS", 1000004, value, entry+" ("+tag+")");
        //          else if(Utils::iequals(entry, "~u_3")) SLHAea_add(out, "MASS", 1000006, value, entry+" ("+tag+")");
        //          else if(Utils::iequals(entry, "~u_4")) SLHAea_add(out, "MASS", 2000002, value, entry+" ("+tag+")");
        //          else if(Utils::iequals(entry, "~u_5")) SLHAea_add(out, "MASS", 2000004, value, entry+" ("+tag+")");
        //          else if(Utils::iequals(entry, "~u_6")) SLHAea_add(out, "MASS", 2000006, value, entry+" ("+tag+")");
        //          else if(Utils::iequals(entry, "~e-_1")) SLHAea_add(out, "MASS", 1000011, value, entry+" ("+tag+")");
        //          else if(Utils::iequals(entry, "~e-_2")) SLHAea_add(out, "MASS", 1000013, value, entry+" ("+tag+")");
        //          else if(Utils::iequals(entry, "~e-_3")) SLHAea_add(out, "MASS", 1000015, value, entry+" ("+tag+")");
        //          else if(Utils::iequals(entry, "~e-_4")) SLHAea_add(out, "MASS", 2000011, value, entry+" ("+tag+")");
        //          else if(Utils::iequals(entry, "~e-_5")) SLHAea_add(out, "MASS", 2000013, value, entry+" ("+tag+")");
        //          else if(Utils::iequals(entry, "~e-_6")) SLHAea_add(out, "MASS", 2000015, value, entry+" ("+tag+")");
        //          else if(Utils::iequals(entry, "~nu_1")) SLHAea_add(out, "MASS", 1000012, value, entry+" ("+tag+")");
        //          else if(Utils::iequals(entry, "~nu_2")) SLHAea_add(out, "MASS", 1000014, value, entry+" ("+tag+")");
        //          else if(Utils::iequals(entry, "~nu_3")) SLHAea_add(out, "MASS", 1000016, value, entry+" ("+tag+")");
        //          else
        //          {
        //             std::ostringstream err;
        //             err << "Error! HDF5Reader encountered an error while attempting to read a spectrum of type '"<<spec_type<<"' from the HDF5 file:group "<<file<<":"<<group<<"' (while calling 'retrieve'). An unrecognised entry with tag 'Pole_Mass' was encountered (label was "<<entry<<")."<<std::endl<<"The full HDF5 entry label was: "<<*it;
        //             printer_error().raise(LOCAL_INFO,err.str());
        //          }
        //       }
        //       else if(Utils::iequals(tag,"mass1"))
        //       {
        //          SLHAea_check_block(out, "MSOFT");
        //          SLHAea_check_block(out, "HMIX");
        //          if     (Utils::iequals(entry, "M1")) SLHAea_add(out, "MSOFT", 1, value, entry+" ("+tag+")");
        //          else if(Utils::iequals(entry, "M2")) SLHAea_add(out, "MSOFT", 2, value, entry+" ("+tag+")");
        //          else if(Utils::iequals(entry, "M3")) SLHAea_add(out, "MSOFT", 3, value, entry+" ("+tag+")");
        //          else if(Utils::iequals(entry, "Mu")) SLHAea_add(out, "HMIX", 1, value, entry+" ("+tag+")");
        //          else if(Utils::iequals(entry, "d_1")) ;
        //          else if(Utils::iequals(entry, "d_2")) ;
        //          else if(Utils::iequals(entry, "d_3")) ;
        //          else if(Utils::iequals(entry, "u_1")) ;
        //          else if(Utils::iequals(entry, "u_2")) ;
        //          else if(Utils::iequals(entry, "u_3")) ;
        //          else if(Utils::iequals(entry, "e-_1")) ;
        //          else if(Utils::iequals(entry, "e-_2")) ;
        //          else if(Utils::iequals(entry, "e-_3")) ;
        //          else if(Utils::iequals(entry, "g")) ;
        //          else if(Utils::iequals(entry, "gamma")) ;
        //          //...etc. TODO: incomplete!
        //       } 
        //     }
        //     else
        //     {
        //       // If one entry is 'invalid' then we cannot reconstruct
        //       // the Spectrum object, so we mark the whole thing invalid.
        //       is_valid = false;
        //     }
        //   }
        // }

        // if(not found_at_least_one)
        // {
        //   // Didn't find any matches!
        //    std::ostringstream err;
        //    err << "Error! HDF5Reader failed to find any Spectrum entries matching the spectrum type '"<<spec_type<<"' in the HDF5 file:group "<<file<<":"<<group<<"' (while calling 'retrieve'). Please check that spectrum type and input file/group are correct.";
        //    printer_error().raise(LOCAL_INFO,err.str());
        // }
        // /// done!
        
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
