#include <iostream>
#include <set>
#include <cstring>

#include "feynrules.hpp"

void FeynRules::load_feynrules()
{

    std::cout << "Loading FeynRules... ";
    
    std::string input;
    // TODO -- somehow the user will need to specify this.
    input+= "$FeynRulesPath = SetDirectory[\"~/.Mathematica/Applications/feynrules-current\"]";
    
    send_to_math(input);
    
    const char* out;
    if (!WSGetString((WSLINK)pHandle, &out))
    {
        std::cout << "Error loading FeynRules. Please check your $FeynRulesPath." << std::endl;
        return;
    }
    else
    {
        std::cout << "FeynRules loaded from " << out << std::endl;
    }
    
    input+= "<<FeynRules`";
    send_to_math(input);
    
}

void FeynRules::load_model(std::string name)
{

    std::cout << "Loading model " + name + " in FeynRules... ";

    // LoadModel command.
    std::string command = "LoadModel[\"Models/" + name + "/" + name + ".fr\"]";
    send_to_math(command);
    
    // Some sort of check here ?
    std::cout << "Model " + name + " loaded successfully." << std::endl;   
     
}

void FeynRules::load_restriction(std::string name)
{

    std::cout << "Loading restriction " + name + "... ";

    // LoadModel command.
    std::string command = "LoadRestriction[\"Models/" + name + ".rst\"]";
    send_to_math(command);
    
    // Some sort of check here ?
    std::cout << "Restriction " + name + " loaded successfully." << std::endl;    
    
}

// Check a model is Hermitian (it should be...)
void FeynRules::check_herm()
{

    std::cout << "Checking the model is Hermitian... ";
    
    // Name of Lagrangian "LTotal" assumed -- probably should be a user input...?
    std::string command = "ch = CheckHermiticity[LTotal]";
    send_to_math(command);
    
    command = "Length[ch]";
    send_to_math(command);
    
    int lench;

    if (!WSGetInteger((WSLINK)pHandle, &lench))
    {
        std::cout << "Error getting 'Length[ch]' from WSTP." << std::endl;
        return;
    }
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
        command = "Length[pl[[" + std::to_string(i+1) + ",2]]]";
        send_to_math(command);
        if (!WSGetInteger((WSLINK)pHandle, &numelements))
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
            int spinX2 = 0; // Needs to be initialised to suppress compiler warnings.
            int pdg;
            bool SM;

            // Assume a particle isn't SC unless we spot it.
            bool self_conjugate = false; 
        
            // Firstly, find the name
            command = "pl[[" + std::to_string(i+1) + ",2," + std::to_string(j+1) + ",1]]";
            send_to_math(command);
            
            if (!WSGetString((WSLINK)pHandle, &name))
            {
                std::cout << "Error getting particle name from WSTP." << std::endl;
                return;
            }
            
            // Then, the antiparticle name
            command = "pl[[" + std::to_string(i+1) + ",2," + std::to_string(j+1) + ",2]]";
            send_to_math(command);
            
            if (!WSGetString((WSLINK)pHandle, &antiname))
            {
                std::cout << "Error getting antiparticle name from WSTP." << std::endl;
                return;
            }
            else if (not (std::strcmp(name, antiname)))
            {
                self_conjugate = true;
            }
            
            // Next, find out the spin.
            command = "pl[[" + std::to_string(i+1) + ",2," + std::to_string(j+1) + ",3]]";
            send_to_math(command);
            
            if (!WSGetString((WSLINK)pHandle, &spin))
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
            if (!WSGetInteger((WSLINK)pHandle, &pdg))
            {
                std::cout << "Error getting pdg code from WSTP." << std::endl;
                return;
            }
            
            // "Full name"
            command = "pl[[" + std::to_string(i+1) + ",2," + std::to_string(j+1) + ",10]]";
            send_to_math(command);
            
            if (!WSGetString((WSLINK)pHandle, &fullname))
            {
                std::cout << "Error getting fullname from WSTP." << std::endl;
                return;
            }
            
            // Check it doesn't get eaten - this isn't a physical particle.
            command = "pl[[" + std::to_string(i+1) + ",2," + std::to_string(j+1) + ",13]]";
            send_to_math(command);
            
            if (!WSGetString((WSLINK)pHandle, &eaten))
            {
                std::cout << "Error getting Goldstone info from WSTP." << std::endl;
                return;
            }
            if (std::strcmp(eaten, "NoGS"))
            {
                continue;
            }
            
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
            Particle particle(pdg, std::string(name), spinX2, std::string(fullname), SM);
            partlist.push_back(particle);
            
            // Also add the antiparticle if it is distinct.
            if (not self_conjugate)
            {
                Particle antiparticle((-1)*pdg, std::string(antiname), spinX2, std::string(fullname), SM);
                partlist.push_back(antiparticle);
            }
            
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

    if (!WSGetInteger((WSLINK)pHandle, &lenepl))
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
        
        if (!WSGetString((WSLINK)pHandle, &block))
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
        if (!WSGetInteger((WSLINK)pHandle, &numelements))
        {
            std::cout << "Error getting number of elements from WSTP." << std::endl;
            return;
        }
        
        for (int j=0; j<numelements; j++) 
        {
            
            const char* paramname;
        
            command = "epl[[" + std::to_string(i+1) + ",2," + std::to_string(j+1) + ",2,1]]";
            send_to_math(command);
            
            if (!WSGetString((WSLINK)pHandle, &paramname))
            {
                std::cout << "Error getting paramname from WSTP." << std::endl;
                return;
            }
            std::cout << "paramname = " << paramname << std::endl;
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
    // Name of Lagrangian "LTotal" assumed -- probably should be a user input...?
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
void FeynRules::write_ch_output(Options opts)
{
    std::cout << "Writing CalcHEP output." << std::endl;
    
    // CalcHEP is faster in Feynman Gauge.
    std::string gauge = "feynman";
    set_gauge(gauge);
    
    // Write output.
    std::string command = "WriteCHOutput[LTotal, CHAutoWidths -> False];";
    send_to_math(command);
    
    /* Checks -- including:
         * Folder where files are written to.
         * Writing of all files.
         * "Done" message.      
    */
    
    std::cout << "CalcHEP files written." << std::endl;
}

// Performs all FeynRules output.
void all_feynrules(Options opts, std::vector<Particle> &partlist, std::vector<Parameter> &paramlist)
{
    std::cout << "Calling FeynRules with model " << opts.model << std::endl;

    // Declare FeynRules model class
    FeynRules model;
 
    // Open link to Mathematica
    model.create_wstp_link();
    
    // Load FeynRules
    model.load_feynrules();
    
    // Set the FeynRules model and load it up
    model.set_name(opts.model);
    model.load_model(opts.model);
    
    // Load restrictions - if there are any.
    if (not opts.restriction.empty())
    {
        model.load_restriction(opts.restriction);
    }
    
    // Diagnositics -- check it is hermitian
    model.check_herm();
    
    // Get all of the particles
    model.get_partlist(partlist);
    
    // And all parameters
    model.get_paramlist(paramlist);
    
    // Write CalcHEP output
    model.write_ch_output(opts);
    
    // All done. Close the Mathematica link.
    model.close_wstp_link();
    
    return;
}
