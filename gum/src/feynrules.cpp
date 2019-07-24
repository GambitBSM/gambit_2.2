//   GUM: GAMBIT Universal Models
//   **********************************
///  \file
///
///  Definitions of Feynrules class
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

#include "feynrules.hpp"

namespace GUM
{

  // Feynrules constructor, loads Feynrules and model
  FeynRules::FeynRules(std::string model, std::string base_model) : Math_Package(model)
  {
    // Math_Package constructor already creates the WSTP link
    
    try
    {
      // Load Feynrules
      load_feynrules();

      // Load model
      load_model(model, base_model);
    }
    catch(...) { throw; }
  }

  // Load Feynrules
  void FeynRules::load_feynrules()
  {

    std::cout << "Loading FeynRules... ";

    std::string input;
    input+= "$FeynRulesPath = SetDirectory[\"" + std::string(FEYNRULES_PATH) + "\"]";

    send_to_math(input);

    const char* out;
    if (!WSGetString(link, &out))
    {
        std::cerr << "Error loading FeynRules. Please check that FeynRules actually lives"
                  << "\nwhere CMake put it, in:\n"
                  << "  " + std::string(FEYNRULES_PATH)
                  << "\nPlease try rebuilding." << std::endl;
        return;
    }
    else
    {
        std::cout << "FeynRules loaded from " << out << "." << std::endl;
    }

    input+= "<<FeynRules`";
    send_to_math(input);

  }

  void FeynRules::load_model(std::string model, std::string base_model)
  {

    std::cout << "Loading model " + model + " in FeynRules";
    if (not base_model.empty()) { std::cout << ", piggybacking off of " << base_model; }
    std::cout << "... " << std::endl;

    // LoadModel command.

    if (base_model.empty()) 
    {
        std::string command = "LoadModel[\"Models/" + model + "/" + model + ".fr\"]";
        send_to_math(command);
    }
    else
    {
        std::string command = "LoadModel[\"Models/" + base_model + "/" + base_model + ".fr\",\"Models/" + model + "/" + model + ".fr\"]";
        send_to_math(command);
    }

    // Check the model has been loaded by querying the model name. If it has changed from the default then we're set.
    // TODO: need to check for duplicate definitions of gauge groups, field contents etc - this makes gum freeze 
    std::string modelname = get_modelname();

    std::string tomatch = "Models" + model + model;
    if (modelname == tomatch)
        throw std::runtime_error("FeynRules Error: Could not load model " + model + ". Please check your FeynRules file.");
    else if ((modelname == "Standard Model") and (base_model.empty()))
        throw std::runtime_error("FeynRules Error: GUM is for BSM physics, yo!");
    else if ((modelname == "Standard Model") and (not base_model.empty()))
        throw std::runtime_error("FeynRules Error: GUM tried to import something on top of the Standard Model, but your model file did not import properly. Please check it.");

    // All good. else {
    std::cout << "Model " + model + " loaded successfully, with model name " << modelname << "." << std::endl;
    //}
  }

  // The model may have a different "internal" name than what's on the package.
  // Need this info for output files, etc.
  std::string FeynRules::get_modelname()
  {
    std::string command = "M$ModelName";
    send_to_math(command);

    const char* out;
    if (!WSGetString(link, &out))
    {
        throw std::runtime_error("FeyRules Error: Error getting Model name.");
    }

    return std::string(out);
  }

  void FeynRules::load_restriction(std::string model, std::string rst)
  {
   
    std::cout << "Loading restriction " + rst + "... ";

    // LoadModel command.
    std::string command = "LoadRestriction[\"Models/" + model + "/" + rst + ".rst\"]";
    std::cout << command << std::endl;
    send_to_math(command);

    // Some sort of check here?
    // TODO: TG: This checks are really clanky
    command = "Length[M$Restrictions]";
    send_to_math(command);

    int out;
    if (!WSGetInteger(link, &out))
        throw std::runtime_error("FeynRules Error: Error loading restriction.");

    if (out == 0)
        throw std::runtime_error("FeynRules Error: No restrictions loaded. Please check your .gum file and \nthat your .rst files are in the correct place.");

    std::cout << "Restriction " + rst + " loaded successfully." << std::endl;
  }

  // Check a model is Hermitian (it should be...)
  void FeynRules::check_herm(std::string LTot)
  {

    std::cout << "Checking the model is Hermitian... ";

    std::cout << "The total Lagrangian is given by : " << LTot << std::endl;

    std::string command = "ch = CheckHermiticity[" + LTot + "]";
    send_to_math(command);

    command = "Length[ch]";
    send_to_math(command);

    int lench;

    if (!WSGetInteger(link, &lench))
        throw std::runtime_error("Error getting 'Length[ch]' from WSTP.");
    else
    {
        if (lench == 0)
        {
            std::cout << "Your Lagrangian is Hermitian." << std::endl;
        }
        else
        {
            std::cout << "Your Lagrangian is not Hermitian." << std::endl;
            std::cout << "FeynRules found " + std::to_string(lench) + " vertices in L-HC[L]." << std::endl;
            std::cout << "Please check your FeynRules file." << std::endl;
            return;
        }

    }

  }

  // Gets the 'PartList' -- this contains all particles in the model.
  void FeynRules::get_partlist(std::vector<Particle> &partlist)
  {

    std::cout << "Extracting particles from FeynRules model." << std::endl;

    // Command to get a list of lists of lists of lists(!) with particle info.
    std::string command = "pl = PartList;";
    send_to_math(command);

    // Find out how many particles we have to get.
    command = "Length[pl]";
    send_to_math(command);

    int lenpl;

    if (!WSGetInteger(link, &lenpl))
    {
        std::cout << "Error getting 'Length[PartList]' from WSTP." << std::endl;
        return;
    }

    std::cout << "Found " << lenpl << " particle sets." << std::endl;

    // Get to parsing this monster.
    for (int i=0; i<lenpl; i++)
    {

        // First things first, check to see if we are dealing with multiplets.
        // e.g., l = (e, mu, tau).
        int numelements;
        command = "Length[pl[[" + std::to_string(i+1) + ",2]]]";
        send_to_math(command);
        if (!WSGetInteger(link, &numelements))
        {
            std::cout << "Error getting number of elements from WSTP." << std::endl;
            return;
        }

        for (int j=0; j<numelements; j++)
        {

            // Initialise all properties we wish to find out about a particle.
            const char* name;
            const char* antiname;
            const char* spin;
            const char* fullname;
            const char* eaten;
            const char* mass;
            // Needs to initialise these to suppress compiler warnings.
            int color = 1;
            int chargeX3 = 0;
            int spinX2 = 0;   
            int pdg;
            bool SM;

            // Firstly, find the name
            command = "pl[[" + std::to_string(i+1) + ",2," + std::to_string(j+1) + ",1]]";
            send_to_math(command);

            if (!WSGetString(link, &name))
            {
                std::cout << "Error getting particle name from WSTP." << std::endl;
                return;
            }

            // Then, the antiparticle name
            command = "pl[[" + std::to_string(i+1) + ",2," + std::to_string(j+1) + ",2]]";
            send_to_math(command);

            if (!WSGetString(link, &antiname))
            {
                std::cout << "Error getting antiparticle name from WSTP." << std::endl;
                return;
            }

            // Next, find out the spin.
            command = "pl[[" + std::to_string(i+1) + ",2," + std::to_string(j+1) + ",3]]";
            send_to_math(command);

            if (!WSGetString(link, &spin))
            {
                std::cout << "Error getting spin from WSTP." << std::endl;
                return;
            }
            if (not std::strcmp(spin, "V"))
            {
                spinX2 = 2;
            }
            else if (not std::strcmp(spin, "F"))
            {
                spinX2 = 1;
            }
            else if (not std::strcmp(spin, "S"))
            {
                spinX2 = 0;
            }
            // Don't try and parse ghosts; it'll give us an error. We also
            // don't care about them from a pheno (GUM) standpoint.
            else if (not std::strcmp(spin, "U"))
            {
                continue;
            }

            // PDG code.
            command = "pl[[" + std::to_string(i+1) + ",2," + std::to_string(j+1) + ",9]]";
            send_to_math(command);
            if (!WSGetInteger(link, &pdg))
            {
                std::cout << "Error getting pdg code from WSTP." << std::endl;
                return;
            }

            // "Full name"
            command = "pl[[" + std::to_string(i+1) + ",2," + std::to_string(j+1) + ",10]]";
            send_to_math(command);

            if (!WSGetString(link, &fullname))
            {
                std::cout << "Error getting fullname from WSTP." << std::endl;
                return;
            }

            // Name of the mass parameter.
            command = "pl[[" + std::to_string(i+1) + ",2," + std::to_string(j+1) + ",5]]";
            send_to_math(command);

            if (!WSGetString(link, &mass))
            {
                std::cout << "Error getting mass from WSTP." << std::endl;
                return;
            }

            // Check it doesn't get eaten - this isn't a physical particle.
            command = "pl[[" + std::to_string(i+1) + ",2," + std::to_string(j+1) + ",13]]";
            send_to_math(command);

            if (!WSGetString(link, &eaten))
            {
                std::cout << "Error getting Goldstone info from WSTP." << std::endl;
                return;
            }
            if (std::strcmp(eaten, "NoGS"))
            {
                continue;
            }

            // Color representation
            command = "pl[[" + std::to_string(i+1) + ",2," + std::to_string(j+1) + ",7]]";
            send_to_math(command);

            const char* colorstr;

            if (!WSGetString(link, &colorstr))
            {
                std::cout << "Error getting Color info from WSTP." << std::endl;
                return;
            }
            if (!std::strcmp(colorstr, "S"))
            {
                color = 1;
            }
            else if (!std::strcmp(colorstr, "T"))
            {
                color = 3;
            }
            else if (!std::strcmp(colorstr, "O"))
            {
                color = 8;
            }

            // Charge

            // Firsly, need to get the name of the Class in FeynRules (in case there's
            // multiple entries -- they will all have the same QNUMBERS.)
            const char* classname;
            command = "pl[[" + std::to_string(i+1) + ",1,2]]";
            send_to_math(command);

            if (!WSGetString(link, &classname))
            {
                std::cout << "Error getting classname from WSTP." << std::endl;
                return;
            }
            std::string str(classname);
            
            // Now see if there is a charge entry
            float charge = 0.;

            command = "N[Q[" + str + "]]";
            send_to_math(command);

            if (!WSGetFloat(link, &charge))
            {
                std::cout << "Error getting charge info from WSTP." << std::endl;
                return;
            }

            chargeX3 = int(round(charge * 3));

            // If we've got this far, our particle is physical.
            // One last thing... check if it's a SM particle or not.
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
            Particle particle(pdg, std::string(name), spinX2, chargeX3, color, std::string(fullname), SM, mass, std::string(antiname));
            partlist.push_back(particle);
        }

     }

  }

  // Gets the 'ParamList' -- this contains all particles in the model.
  void FeynRules::get_paramlist(std::vector<Parameter> &paramlist)
  {
    std::cout << "Extracting external parameters from FeynRules model." << std::endl;

    std::string command = "epl = EParamList;";
    send_to_math(command);

    // Find out how many parameters we have to get.
    command = "Length[epl]";
    send_to_math(command);

    int lenepl;

    if (!WSGetInteger(link, &lenepl))
    {
        std::cout << "Error getting 'Length[EParamList]' from WSTP." << std::endl;
        return;
    }

    std::cout << "Found " << lenepl << " parameter blocks." << std::endl;

    const char* block;

    for (int i=0; i<lenepl; i++)
    {
        command = "epl[[" + std::to_string(i+1) + ",1]]";
        send_to_math(command);

        if (!WSGetString(link, &block))
        {
            std::cout << "Error getting block info from WSTP. "
                      << "Please check every parameter has a BlockName in your .fr file."
                      << std::endl;
            return;
        }

        // First things first, check how many elements are in each block.
        // e.g., SMINPUTS = (aEWM1, Gf, alphaS).
        int numelements;
        command = "Length[epl[[" + std::to_string(i+1) + ",2]]]";
        send_to_math(command);
        if (!WSGetInteger(link, &numelements))
        {
            std::cout << "Error getting number of elements from WSTP." << std::endl;
            return;
        }

        for (int j=0; j<numelements; j++)
        {

            const char* paramname;

            command = "epl[[" + std::to_string(i+1) + ",2," + std::to_string(j+1) + ",2,1]]";
            send_to_math(command);

            if (!WSGetString(link, &paramname))
            {
                std::cout << "Error getting paramname from WSTP." << std::endl;
                return;
            }

            Parameter parameter(paramname, block, int(j+1));
            paramlist.push_back(parameter);
        }
    }
  }

  // Set the gauge -- accepts either 'feynman' or 'unitary'
  void FeynRules::set_gauge(std::string gauge)
  {
    if ((gauge != "feynman") && (gauge != "unitary"))
    {
        std::cerr << "set_gauge() called with an incorrect option, "
                  << gauge << ".\n Allowed options: 'feynman' or 'unitary'."
                  << std::endl;
        return;
    }
    std::string command;
    if (gauge == "feynman")
    {
        command = "FeynmanGauge = True;";
        std::cout << "Setting Feynman Gauge." << std::endl;
    }
    else if (gauge == "unitary")
    {
        command = "FeynmanGauge = False;";
        std::cout << "Setting Unitary Gauge." << std::endl;
    }
    send_to_math(command);
  }

  // Write CalcHEP output.
  void FeynRules::write_ch_output(std::string LTot)
  {
    std::cout << "Writing CalcHEP output." << std::endl;

    // CalcHEP is faster in Feynman Gauge.
    std::string gauge = "feynman";
    set_gauge(gauge);

    // Write output.
    std::string command = "WriteCHOutput[" + LTot + ", CHAutoWidths -> False];";
    send_to_math(command);

    std::cout << "CalcHEP files written." << std::endl;
  }

  // Write MadGraph output.
  void FeynRules::write_mg_output(std::string LTot)
  {
    std::cout << "Writing MadGraph (UFO) output for Pythia/MadDM." << std::endl;

    // Write output.
    std::string command = "WriteUFO[" + LTot + "];";
    send_to_math(command);

    std::cout << "MadGraph files written." << std::endl;
  }

  // Performs all FeynRules output.
  void all_feynrules(Options opts, std::vector<Particle> &partlist, std::vector<Parameter> &paramlist, Outputs &outputs, std::vector<std::string> &backends)
  {
    try
    {
      std::cout << "Calling FeynRules with model " << opts.model() << "..." << std::endl;

      // Create FeynRules object, open link to Mathematica, load Feynrules and the model
      FeynRules model(opts.model(), opts.base_model());

      // Load restrictions - if there are any.
      if (not opts.restriction().empty())
      {
        model.load_restriction(opts.model(), opts.restriction());
      }

      // Diagnositics -- check it is hermitian
      model.check_herm(opts.lagrangian());

      // Get all of the particles
      model.get_partlist(partlist);

      // And all parameters
      model.get_paramlist(paramlist);

      // Get the "actual" model name
      std::string fr_model_name = model.get_modelname();

      // Output directory
      std::string output = std::string(FEYNRULES_PATH) + "/" + fr_model_name;

      /// Write CalcHEP output
      if (std::find(backends.begin(), backends.end(), "calchep") != backends.end() )
        model.write_ch_output(opts.lagrangian());

      // Location of CalcHEP files
      std::string chdir = output + "-CH";
      std::replace(chdir.begin(), chdir.end(), ' ', '-');
      outputs.set_ch(chdir);

      /// Same for MadGraph->Pythia
      if (std::find(backends.begin(), backends.end(), "pythia") != backends.end() )
        model.write_mg_output(opts.lagrangian());

      // Location of MadGraph (UFO) files
      std::string mgdir = output + "_UFO";
      std::replace(mgdir.begin(), mgdir.end(), ' ', '_');
      outputs.set_mg(mgdir);

    }
    catch(std::exception &e)
    {
     std::cerr << e.what() << std::endl;
    }
  }

} // namespace GUM

// Now all the grizzly stuff, so Python can call C++ (which can call Mathematica...)
BOOST_PYTHON_MODULE(libfr)
{
  using namespace boost::python;

  class_<Particle>("FRParticle", init<int, std::string, int, int, int, std::string, bool, std::string, std::string>())
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

  class_<Parameter>("FRParameter", init<std::string, std::string, int>())
    .def("name",  &Parameter::name)
    .def("block", &Parameter::block)
    .def("index", &Parameter::index)
    ;

  class_<Options>("FROptions", init<std::string, std::string, std::string, std::string, std::string>())
    .def("package",     &Options::package)
    .def("model",       &Options::model)
    .def("base_model",  &Options::base_model)
    .def("restriction", &Options::restriction)
    .def("lagrangian",  &Options::lagrangian)
    ;

  class_<Outputs>("FROutputs", init<>())
    .def("get_ch",  &Outputs::get_ch)
    .def("get_mg",  &Outputs::get_mg)
    ;

  class_< std::vector<Particle> >("FRVectorOfParticles")
    .def(vector_indexing_suite< std::vector<Particle> >() )
    ;

  class_< std::vector<Parameter> >("FRVectorOfParameters")
    .def(vector_indexing_suite< std::vector<Parameter> >() )
    ;

  class_< std::vector<std::string> >("FRBackends")
    .def(vector_indexing_suite< std::vector<std::string> >() )
    ;


  def("all_feynrules", GUM::all_feynrules);

}
