#include <iostream>
#include <set>
#include <algorithm>
#include <cstring>

#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

#include "sarah.hpp"

void SARAH::load_sarah()
{

    std::cout << "Loading SARAH... ";
    
    std::string input;       
    input+= "SetDirectory[\"" + std::string(SARAH_PATH) + "\"]";        
    send_to_math(input);
    
    const char* out;
    if (!WSGetString((WSLINK)pHandle, &out))
    {
        std::cerr << "Error loading SARAH. Please check that SARAH actually lives" 
                  << "\nwhere CMake put it, in:\n" 
                  << "  " + std::string(SARAH_PATH)
                  << "\nPlease try rebuilding." << std::endl;
        return;
    }
    else
    {
        std::cout << "SARAH loaded from " << out << "." << std::endl;
    }
    
    input+= "<<SARAH`";
    send_to_math(input);
    
}

bool SARAH::load_model(std::string model)
{

    std::cout << "Loading model " + model + " in SARAH... " << std::endl;

    // Load it up.
    std::string command = "Start[\"" + model + "\"];";
    send_to_math(command);
    
    // Check the model has been loaded by querying the model name. If it has changed from the default then we're set.
    std::string modelname;
    get_modelname(modelname);

    // ...Assuming someone hasn't set the model name to 'ModelName' which would be unbelievably annoying and vastly silly.
    if (modelname == "ModelName")
    {
        std::cerr << std::endl << "ERROR! Could not load model " << model << ". Please check your SARAH file." << std::endl << std::endl;
        return false;
    }
    
    // All good.
    std::cout << "Model " + model + " loaded successfully, with model name " << modelname << "." << std::endl;  
    return true;
}

// The model may have a different "internal" name than what's on the package.
// Need this info for output files, etc.
void SARAH::get_modelname(std::string &modelname)
{

    std::string command = "ModelName";
    send_to_math(command);
    
    const char* out;
    if (!WSGetString((WSLINK)pHandle, &out))
    {
        std::cerr << "Error getting Model name." << std::endl;
        return;
    }
    
    modelname = std::string(out);
}

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
    
    if (!WSGetInteger((WSLINK)pHandle, &lenpl))
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
        command = "Length[getPDG[pl[[" + std::to_string(i+1) + ", 1]]]]";
        send_to_math(command);
        if (!WSGetInteger((WSLINK)pHandle, &numelements))
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
            const char* fullname;
            std::string mass;
            int spinX2 = 0; // Needs to be initialised to suppress compiler warnings.
            int pdg;
            int num;
            bool SM;
            bool capitalise = false;
            
            // Assume a particle is SC unless we spot it.
            bool self_conjugate = true; 
            
            command = "Part[getPDG[pl[[" + std::to_string(i+1) + ", 1]]], " + std::to_string(j+1) + "]";
            send_to_math(command);
            
            if (!WSGetInteger((WSLINK)pHandle, &pdg))
            {
                std::cout << "Error getting PDG code from WSTP." << std::endl;
                return;
            }
            
            // If it's got a PDG of 0 it's not a physical particle. Don't care about it.
            if (pdg == 0) { continue; }
            
            command = "Length[getOutputName[pl[[" + std::to_string(i+1) + ", 1]]]]";
            send_to_math(command);
            
            if (!WSGetInteger((WSLINK)pHandle, &num))
            {
                std::cout << "Error getting length of OutputNames." << std::endl;
                return;
            }
            
            if (num == 2) 
            {
                self_conjugate = false;        
                command = "Part[getOutputName[pl[[" + std::to_string(i+1) + ", 1]]], 1]";
                send_to_math(command);
                
                if (!WSGetString((WSLINK)pHandle, &name))
                {
                    std::cerr << "Error getting particle name." << std::endl;
                    return;
                }
                command = "Part[getOutputName[pl[[" + std::to_string(i+1) + ", 1]]], 2]";
                send_to_math(command);
                
                if (!WSGetString((WSLINK)pHandle, &antiname))
                {
                    std::cerr << "Error getting particle antiname." << std::endl;
                    return;
                }
            }
            else if (num == 0)
            {
                command = "getOutputName[pl[[" + std::to_string(i+1) + ", 1]]]";
                send_to_math(command);
                
                if (!WSGetString((WSLINK)pHandle, &name))
                {
                    std::cerr << "Error getting particle name." << std::endl;
                    return;
                }
                                
                // Probe to see if it is self-conjugate
                command = "TrueQ[pl[[" + std::to_string(i+1) + ", 1]] == conj[pl[[" + std::to_string(i+1) + ", 1]]]]";
                send_to_math(command);
                
                const char* is_sc;
                if (!WSGetString((WSLINK)pHandle, &is_sc))
                {
                    std::cerr << "Error getting name." << std::endl;
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
                std::cerr << "More than 2 particles here; what weird symmetries have you got???" << std::endl;
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
            }
            
            //std::cout << "PDG code = " << pdg << ", self conj? -- " << self_conjugate << ", outputname = " << outputname;
            //if (not self_conjugate) { std::cout << ", antioutputname = " << antioutputname << std::endl; }
            //else { std::cout << std::endl; }
            
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
            Particle particle(pdg, std::string(name), spinX2, std::string(outputname), SM, mass, self_conjugate);
            partlist.push_back(particle);
                        
            // Also add the antiparticle if it is distinct.
            if (not self_conjugate)
            {
                Particle antiparticle((-1)*pdg, std::string(antiname), spinX2, std::string(antioutputname), SM, mass, self_conjugate);
                partlist.push_back(antiparticle);
            }
                                    
        }
        
    }
             
}

void SARAH::get_paramlist(std::vector<Parameter> &paramlist)
{

    std::cout << "Extracting parameters from SARAH model." << std::endl;
    
    // Get list of parameters
    std::string command = "pL = parameters;";
    send_to_math(command);

    // Find out how many parameters we have to get.
    command = "Length[pL]";
    send_to_math(command);
    
    int lenpl;
    
    if (!WSGetInteger((WSLINK)pHandle, &lenpl))
    {
        std::cout << "Error getting 'Length[ParamList]' from WSTP." << std::endl;
        return;
    }
    
    std::cout << "Found " << lenpl << " parameter sets." << std::endl;

}

// Write CalcHEP output. 
void SARAH::write_ch_output()
{
    std::cout << "Writing CalcHEP output." << std::endl;
        
    // Write output.
    std::string command = "MakeCHep[SLHAinput -> False, UseRunningCoupling -> True, CalculateMasses -> False];";
    send_to_math(command);
        
    std::cout << "CalcHEP files written." << std::endl;
}

void all_sarah(Options opts, std::vector<Particle> &partlist, std::vector<Parameter> &paramlist, Outputs &outputs)
{

    std::cout << "Calling SARAH with model " << opts.model() << "..." << std::endl;
    
    // Initialise SARAH object
    SARAH model;
    
    // Open the WSTP link to Mathematica
    model.create_wstp_link();
    
    // Load SARAH.
    model.load_sarah();
    
    // Set the name for the model class.
    model.set_name(opts.model());
    
    // Attempt to load the SARAH model
    bool out = model.load_model(opts.model());
    
    if (not out)
    {
        return;
    }
    
    // Get all of the particles
    model.get_partlist(partlist);
    
    // And all parameters
    model.get_paramlist(paramlist);
    
    // Write CalcHEP output
    model.write_ch_output();
    
    // This is where the output CalcHEP files will be written
    std::string chdir = std::string(SARAH_PATH) + "/Output/" + opts.model() + "/EWSB/CHep";
    std::replace(chdir.begin(), chdir.end(), ' ', '-');
    
    outputs.set_ch(chdir);
    
    // All done. Close the Mathematica link.
    model.close_wstp_link();
    
    return;
    
}

// Now all the grizzly stuff, so Python can call C++ (which can call Mathematica...) 
BOOST_PYTHON_MODULE(libsarah)
{
  using namespace boost::python;
    
  class_<Particle>("SARAHParticle", init<int, std::string, int, std::string, bool, std::string, bool>())
    .def("pdg",    &Particle::pdg)
    .def("name",   &Particle::name)
    .def("SM",     &Particle::SM)
    .def("spinX2", &Particle::spinX2)
    .def("mass",   &Particle::mass)
    .def("SC",     &Particle::SC)
    ;
  
  class_<Parameter>("SARAHParameter", init<std::string, std::string>())
    .def("name",  &Parameter::name)
    .def("block", &Parameter::block)
    ;
  
  class_<Options>("SARAHOptions", init<std::string, std::string, std::string>())
    .def("package",     &Options::package)
    .def("model",       &Options::model)
    ;
  
  class_<Outputs>("SARAHOutputs", init<>())
    .def("get_ch",  &Outputs::get_ch)
    ;
  
  class_< std::vector<Particle> >("SARAHVectorOfParticles")
    .def(vector_indexing_suite< std::vector<Particle> >() )
    ;
    
  class_< std::vector<Parameter> >("SARAHVectorOfParameters")
    .def(vector_indexing_suite< std::vector<Parameter> >() )
    ;
    

  def("all_sarah", all_sarah);
  
}
