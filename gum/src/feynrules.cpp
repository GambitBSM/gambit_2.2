#include "feynrules.hpp"
#include <iostream>

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

// Check a model is Hermitian (it should be...)
void FeynRules::check_herm()
{

    std::cout << "Checking the model is Hermitian... ";
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
