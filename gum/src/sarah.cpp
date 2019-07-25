//   GUM: GAMBIT Universal Models
//   **********************************
///  \file
///
///  Definitions of SARAH class
///
///  **********************************
///
///  \author Sanjay Bloor
///          (sanjay.bloor12@imperial.ac.uk)
///  \date 2017, 2018, 2019
///
///  \author Tomas Gonzalo
///          (tomas.gonzalo@monash.edu)
///  \date 2019 July
///
///  ***********************************

#include <set>
#include <algorithm>
#include <cstring>

#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

#include "sarah.hpp"

namespace GUM
{

  // SARAH constructor, loads SARAH and the model
  SARAH::SARAH(std::string model) : Math_Package(model)
  {
    // Math_Package constructor already creates the WSTP link

    try
    {
      // Load SARAH
      load_sarah();

      // Load model
      load_model(model);
    }
    catch(...) { throw; }
    
  }
   
  // Load SARAH
  void SARAH::load_sarah()
  {

    std::cout << "Loading SARAH... ";

    std::string input;
    input+= "SetDirectory[\"" + std::string(SARAH_PATH) + "\"]";
    send_to_math(input);

    const char* out;
    if (!WSGetString(link, &out))
        throw std::runtime_error("SARAH Error: Error loading SARAH. Please check that SARAH actually lives \nwhere CMake put it, in:\n " + std::string(SARAH_PATH) + "\nPlease try rebuilding.");
    else
    {
        std::cout << "SARAH loaded from " << out << "." << std::endl;
    }

    input+= "<<SARAH`";
    send_to_math(input);

  }

  void SARAH::load_model(std::string model)
  {
    try
    {
      std::cout << "Loading model " + model + " in SARAH... " << std::endl;

      // Check if model is in SARAH's or GUM's list of models
      if(!check_model(model))
      {
        throw std::runtime_error("SARAH Error: Could not load model " + model + ". Model is not recognised by SARAH or GUM.");
      }

      // Load it up.
      std::string command = "Start[\"" + model + "\"];";
      send_to_math(command);

      // Check the model has been loaded by querying the model name. If it has changed from the default then we're set.
      std::string modelname = get_modelname();

      // ...Assuming someone hasn't set the model name to 'ModelName' which would be unbelievably annoying and vastly silly.
      if (modelname == "ModelName")
        throw std::runtime_error("SARAH Error: Could not load model " + model + ". Please check your SARAH file.");


      // All good.
      std::cout << "Model " + model + " loaded successfully, with model name " << modelname << "." << std::endl;
  
    } catch(...) { throw; }
  }

  // The model may have a different "internal" name than what's on the package.
  // Need this info for output files, etc.
  std::string SARAH::get_modelname()
  {

    std::string command = "ModelName";
    send_to_math(command);

    const char* out;
    if (!WSGetString(link, &out))
        throw std::runtime_error("SARAH Error: Error getting Model name.");

    return std::string(out);
  }

  // Check if model is SARAH's database or in GUM's
  bool SARAH::check_model(std::string modelname)
  {
    try
    {
      // Check if the model is in the SARAH database
      std::string command = "MemberQ[ShowModels[[1]],\"" + modelname + "\"]";
      send_to_math(command);
      
      // Get the boolean result
      bool is_SARAH_model;
      get_from_math(is_SARAH_model);
      if(is_SARAH_model)
        return true;

      // If not check if it's on the GUM model database
      std::string model_paths = std::string(GUM_DIR) + "/Models";
      // TODO: check if it's in the models dir and if it is move it to the SARAH dir

      return false;
    
    }
    catch(...) { throw; }
  }

  // Get particles list
  void SARAH::get_partlist(std::vector<Particle> &partlist)
  {

    std::cout << "Extracting particles from SARAH model." << std::endl;

    // Command to get a list with (most) particle info.
    std::string command = "pl = ParticleDefinitions[EWSB];";
    send_to_math(command);

    // Find out how many particles we have to get.
    command = "Length[pl]";
    send_to_math(command);

    int lenpl;

    if (!WSGetInteger(link, &lenpl))
        throw std::runtime_error("SARAH Error: Error getting 'Length[PartList]' from WSTP.");

    std::cout << "Found " << lenpl << " particle sets." << std::endl;

    // Get to parsing this monster.
    for (int i=0; i<lenpl; i++)
    {

        // First things first, check to see if we are dealing with multiplets.
        // e.g., l = (e, mu, tau).
        int numelements;
        command = "Length[getPDG[pl[[" + std::to_string(i+1) + ", 1]]]]";
        send_to_math(command);
        if (!WSGetInteger(link, &numelements))
        {
            std::cout << "Error getting number of elements from WSTP." << std::endl;
            return;
        }

        // If there's no associated PDG code.
        if (numelements == 0)
        {
            continue;
        }

        for (int j=0; j<numelements; j++)
        {
            // Initialise all properties we wish to find out about a particle.
            const char* name;
            std::string outputname;
            const char* antiname;
            std::string antioutputname;
            std::string mass;
            int spinX2 = 0; // Needs to be initialised to suppress compiler warnings.
            int chargeX3 = 0;
            int color = 0;
            int pdg;
            int num;
            bool SM;
            bool capitalise = false;

            // Assume a particle is SC unless we spot it.
            bool self_conjugate = true;

            command = "Part[getPDG[pl[[" + std::to_string(i+1) + ", 1]]], " + std::to_string(j+1) + "]";
            send_to_math(command);

            if (!WSGetInteger(link, &pdg))
            {
                std::cout << "Error getting PDG code from WSTP." << std::endl;
                return;
            }

            // If it's got a PDG of 0 it's not a physical particle. Don't care about it.
            if (pdg == 0) { continue; }

            command = "Length[getOutputName[pl[[" + std::to_string(i+1) + ", 1]]]]";
            send_to_math(command);

            if (!WSGetInteger(link, &num))
            {
                std::cout << "Error getting length of OutputNames." << std::endl;
                return;
            }

            if (num == 2)
            {
                self_conjugate = false;
                command = "Part[getOutputName[pl[[" + std::to_string(i+1) + ", 1]]], 1]";
                send_to_math(command);

                if (!WSGetString(link, &name))
                {
                    std::cerr << "Error getting particle name." << std::endl;
                    return;
                }
                command = "Part[getOutputName[pl[[" + std::to_string(i+1) + ", 1]]], 2]";
                send_to_math(command);

                if (!WSGetString(link, &antiname))
                {
                    std::cerr << "Error getting particle antiname." << std::endl;
                    return;
                }
            }
            else if (num == 0)
            {
                command = "getOutputName[pl[[" + std::to_string(i+1) + ", 1]]]";
                send_to_math(command);

                if (!WSGetString(link, &name))
                {
                    std::cerr << "Error getting particle name." << std::endl;
                    return;
                }

                // Probe to see if it is self-conjugate
                command = "TrueQ[pl[[" + std::to_string(i+1) + ", 1]] == conj[pl[[" + std::to_string(i+1) + ", 1]]]]";
                send_to_math(command);

                const char* is_sc;
                if (!WSGetString(link, &is_sc))
                {
                    std::cerr << "Error getting self-conjugate status " 
                              << "for particle " 
                              << std::string(name) << "."
                              << std::endl;
                    return;
                }
                if (strcmp(is_sc, "True"))
                {
                    self_conjugate = false;
                    capitalise = true;
                }
            }
            else
            {
                std::cerr << "More than 2 particles here; "
                          << "what weird symmetries have you got???" 
                          << std::endl;
                return;
            }

            if (numelements > 1)
            {
                outputname = std::string(name) + std::to_string(j+1);
                if (not self_conjugate && capitalise)
                {
                    antioutputname = outputname;
                    if (isupper(antioutputname[0])) { antioutputname = tolower(antioutputname[0]); }
                    else { antioutputname[0] = toupper(antioutputname[0]); }
                }
            }
            else
            {

                outputname = std::string(name);
                if (not self_conjugate)
                {
                    antioutputname = std::string(antiname);
                }
                else
                {
                    antioutputname = std::string(name);
                }
            }

            mass = "M" + outputname;

            std::set<int> SM_pdgs = {1, 2, 3, 4, 5, 6, 11, 12, 13, 14, 15, 16, 21, 22, 23, 24};
            if (SM_pdgs.count(abs(pdg)))
            {
                SM = true;
            }
            else
            {
                SM = false;
            }

            // Add the particle to the list.
            Particle particle(pdg, std::string(outputname), spinX2, chargeX3, color, std::string(outputname), SM, mass, std::string(antioutputname));
            partlist.push_back(particle);

        }

    }

    std::cout << "Done." << std::endl;

  }

  // Get parameters list
  void SARAH::get_paramlist(std::vector<Parameter> &paramlist)
  {

    std::cout << "Extracting parameters from SARAH model." << std::endl;

    // Get list of parameters
    std::string command = "pd = ParameterDefinitions;";
    send_to_math(command);

    // Find out how many parameters we have to get.
    command = "Length[pd]";
    send_to_math(command);

    int lenpl;

    if (!WSGetInteger(link, &lenpl))
    {
        std::cout << "Error getting 'Length[ParamList]' from WSTP." << std::endl;
        return;
    }

    std::cout << "Found " << lenpl << " parameter sets." << std::endl;

    for (int i=0; i<lenpl; i++)
    {
        const char* block;
        const char* paramname;
        int index;

        // Whether or not the parameter is defined by other parameters
        // of the model. If so, don't want it in GAMBIT
        bool externalparam = true;

        // Whether we've found an LH block
        bool LHblock = false;

        command = "pd[[" + std::to_string(i+1) + ",1]]";
        send_to_math(command);

        // Get the parameter name as it is known in SARAH. This
        // might change later.
        if (!WSGetString(link, &paramname))
        {
            std::cout << "Error getting the parameter name "
                      << "at position " << i+1 
                      << "; please check your .m file."
                      << std::endl;
            return;
        }

        command = "Length[pd[[" + std::to_string(i+1) + ",2]]]";
        send_to_math(command);

        int numelements;

        // Find out how many entries are in the description 
        // of the parameter
        if (!WSGetInteger(link, &numelements))
        {
            std::cout << "Error getting the number of elements "
                      << "defining the parameter " << paramname
                      << ". Please check your .m file."
                      << std::endl;
            return;
        }

        // Go through each parameter and extract any useful info
        for (int j=0; j<numelements; j++)
        {

            // Each entry will be some sort of descriptor for
            // a particle. Discard things like LaTeX entries, 
            // descriptions in prose, etc.
            const char* entry;

            // Here are the useful things we want to query:
            // - blockname & index
            // - parameter name
            // - dependences

            // Does the parameter depend on some other combination
            // of parameters?
            command = "DependenceNum /. pd[[" + std::to_string(i+1) + ",2," + std::to_string(j+1) + "]] // ToString";
            send_to_math(command);

            if (!WSGetString(link, &entry))
            {
                std::cout << "Error querying DependenceNum for "
                          << i+1 << ", " << j+1 << " from WSTP "
                          << "for the SARAH parameter " << paramname
                          << std::endl;
                return;
            }
            // If it has a dependence -- not interested. Bin it.
            if (strcmp(entry, "DependenceNum")) 
            { 
                externalparam = false; 
            }

            // Does the parameter depend on some other combination
            // of parameters?
            command = "Head[LesHouches /. pd[[" + std::to_string(i+1) + ",2," + std::to_string(j+1) + "]]]";
            send_to_math(command);

            if (!WSGetString(link, &entry))
            {
                std::cout << "Error querying LesHouches entry "
                          << i+1 << ", " << j+1 << " from WSTP "
                          << "for the SARAH parameter " << paramname
                          << std::endl;
                return;
            }

            // If we have a list, then there's a blockname and an index.
            if (not strcmp(entry, "List")) 
            { 
                // blockname
                command = "Part[LesHouches /. pd[[" + std::to_string(i+1) + ",2," + std::to_string(j+1) + "]],1]";
                send_to_math(command);

                if (!WSGetString(link, &block))
                {
                    std::cout << "Error querying blockname from WSTP "
                              << "for the SARAH parameter " 
                              << paramname << " with indices "
                              << i+1 << "," << j+1
                              << std::endl;
                    return;
                }

                // index
                command = "Part[LesHouches /. pd[[" + std::to_string(i+1) + ",2," + std::to_string(j+1) + "]],2]";
                send_to_math(command);
                if (!WSGetInteger(link, &index))
                {
                    std::cout << "Error querying index from WSTP "
                              << "for the SARAH parameter " 
                              << paramname << " with indices "
                              << i+1 << "," << j+1
                              << std::endl;
                    return;
                }
                LHblock = true;
            }

            // Does the parameter have a different external
            // name than the internal SARAH name? If so, use it.
            command = "OutputName /. pd[[" + std::to_string(i+1) + ",2," + std::to_string(j+1) + "]] // ToString";
            send_to_math(command);

            if (!WSGetString(link, &entry))
            {
                std::cout << "Error querying OutputName for "
                          << i+1 << ", " << j+1 << " from WSTP "
                          << "for the SARAH parameter " << paramname
                          << std::endl;
                return;
            }
            // If it has a different output name to 
            // internal, we'd better use it, seeing as GAMBIT
            // is... output, I guess.
            if (strcmp(entry, "OutputName")) 
            { 
                paramname = entry;
            }
        }

        // If it's a fundamental parameter of our theory, 
        // add it to the list (if it has an LH block!)
        if (externalparam && LHblock)
        {
            Parameter parameter(paramname, block, index);
            paramlist.push_back(parameter);
        }
    }

  }

  // Write CalcHEP output.
  void SARAH::write_ch_output()
  {
    std::cout << "Writing CalcHEP output." << std::endl;

    // Options for the CH output.
    std::string options;
    // This is currently hard-coded, because GAMBIT needs to interface with CalcHEP 
    // in a specific way.
    /* do not read from SLHA file |  let alphaS run  |  all masses computed externally */
    options = "SLHAinput -> False, UseRunningCoupling -> True, CalculateMasses -> False";

    // Write output.
    std::string command = "MakeCHep[" + options + "];";
    send_to_math(command);

    std::cout << "CalcHEP files written." << std::endl;
  }

  // Write MadGraph output.
  void SARAH::write_madgraph_output()
  {
    std::cout << "Writing MadGraph (UFO) output for Pythia/MadDM." << std::endl;

    // Write output.
    std::string command = "MakeUFO[];";
    send_to_math(command);

    std::cout << "MadGraph files written." << std::endl;
  }

  // Write SPheno output.
  void SARAH::write_spheno_output()
  {
      std::cout << "Writing SPheno output." << std::endl;
      
      // Options for SPheno output.
      std::string options;
      // TODO: options:
      // - InputFile (default $MODEL/SPheno.m)
      // - StandardCompiler -> <COMPILER> (default gfortran) // TG: This should be handled by GM cmake system, so no need

      // Write output.
      std::string command = "MakeSPheno[" + options + "];";
      send_to_math(command);

      std::cout << "SPheno files written." << std::endl;
  }

  // Write Vevacious output.
  void SARAH::write_vevacious_output()
  {
      std::cout << "Writing Vevacious output." << std::endl;  

      // Options for Vevacious output.
      std::string options;
      // TODO: options:
      // - ComplexParameters (automatic?)
      // - Scheme (DRbar for SUSY, MSbar for non-SUSY)  

      // Write output.
      std::string command = "MakeVevacious[" + options + "];";
      send_to_math(command);  

      std::cout << "Vevacious files written." << std::endl;
  }

  // Do all operations with SARAH
  void all_sarah(Options opts, std::vector<Particle> &partlist, std::vector<Parameter> &paramlist, Outputs &outputs, std::vector<std::string> &backends)
  {

    try
    {
      std::cout << "Calling SARAH with model " << opts.model() << "..." << std::endl;

      // Create SARAH object, open link to Mathematica, load SARAH and the model
      SARAH model(opts.model());

      // Get all of the particles
      model.get_partlist(partlist);

      // And all parameters
      model.get_paramlist(paramlist);

      // Where the outputs all live
      std::string outputdir = std::string(SARAH_PATH) + "/Output/" + opts.model() + "/EWSB/";

      /// Write CalcHEP output
      if (std::find(backends.begin(), backends.end(), "calchep") != backends.end() )
      {
        model.write_ch_output();

        // Location of CalcHEP files
        std::string chdir = outputdir + "CHep";
        std::replace(chdir.begin(), chdir.end(), ' ', '-');
        outputs.set_ch(chdir);
      }
 
      /// Write MadGraph output
      if (std::find(backends.begin(), backends.end(), "pythia") != backends.end() )
      {
        model.write_madgraph_output();

        // Location of MadGraph (UFO) files
        std::string mgdir = outputdir + "UFO";
        std::replace(mgdir.begin(), mgdir.end(), ' ', '-');
        outputs.set_mg(mgdir);
      }

      /// Write SPheno output
      if (std::find(backends.begin(), backends.end(), "spheno") != backends.end() )
      {
        // TODO: commented for speed, uncomment at the end
        model.write_spheno_output();

        // Location of SPheno files
        std::string sphdir = outputdir + "SPheno";
        std::replace(sphdir.begin(), sphdir.end(), ' ', '-');
        outputs.set_sph(sphdir);
      }

      /// Write Vevacious output
      if (std::find(backends.begin(), backends.end(), "vevacious") != backends.end() )
      {
        model.write_vevacious_output();        // Location of Vevacious (vin) files
        std::string vevdir = outputdir + "Vevacious";
        std::replace(vevdir.begin(), vevdir.end(), ' ', '-');
        outputs.set_vev(vevdir);
      }

    }
    catch(std::exception &e)
    {
      std::cerr << e.what() << std::endl;
    }
  }

   
} // namespace GUM

// Now all the grizzly stuff, so Python can call C++ (which can call Mathematica...)
BOOST_PYTHON_MODULE(libsarah)
{
  using namespace boost::python;

  class_<Particle>("SARAHParticle", init<int, std::string, int, int, int, std::string, bool, std::string, std::string>())
    .def("pdg",      &Particle::pdg)
    .def("name",     &Particle::name)
    .def("SM",       &Particle::SM)
    .def("spinX2",   &Particle::spinX2)
    .def("chargeX3", &Particle::chargeX3)
    .def("color",    &Particle::color)
    .def("mass",     &Particle::mass)
    .def("SC",       &Particle::SC)
    .def("antiname", &Particle::antiname)
    ;

  class_<Parameter>("SARAHParameter", init<std::string, std::string, int>())
    .def("name",  &Parameter::name)
    .def("block", &Parameter::block)
    .def("index", &Parameter::index)
    ;

  class_<Options>("SARAHOptions", init<std::string, std::string, std::string, std::string, std::string>())
    .def("package",     &Options::package)
    .def("model",       &Options::model)
    ;

  class_<Outputs>("SARAHOutputs", init<>())
    .def("get_ch",   &Outputs::get_ch)
    .def("get_mg",   &Outputs::get_mg)
    .def("get_sph",  &Outputs::get_sph)
    .def("get_vev",  &Outputs::get_vev)
    ;

  class_< std::vector<Particle> >("SARAHVectorOfParticles")
    .def(vector_indexing_suite< std::vector<Particle> >() )
    ;

  class_< std::vector<Parameter> >("SARAHVectorOfParameters")
    .def(vector_indexing_suite< std::vector<Parameter> >() )
    ;

  class_< std::vector<std::string> >("SARAHBackends")
    .def(vector_indexing_suite< std::vector<std::string> >() )
    ;

  def("all_sarah", GUM::all_sarah);

}
