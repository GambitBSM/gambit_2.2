#include <iostream>

#include "feynrules.hpp"

void FeynRules::load_feynrules()
{
    std::cout << "Loading FeynRules... ";
    
    std::string input;
    // TODO -- somehow the user will need to specify this.
    input+= "$FeynRulesPath = SetDirectory[\"~/.Mathematica/Applications/feynrules-current\"]";
    
    send_to_math(input);
    
    const char* out;
    if(!WSGetString((WSLINK)pHandle, &out))
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

    if(!WSGetInteger((WSLINK)pHandle, &lench))
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

// Performs all FeynRules output.
void all_feynrules(Options opts)
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
    
    // All done. Close the Mathematica link.
    model.close_wstp_link();
    
    return;
}
