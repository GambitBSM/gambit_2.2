/// Contains all routines for interface between C++ and Mathematica packages.

#include <iostream>

#include "feynrules.hpp"

int main()
{
    std::string frmodel = "SingletDM";
    std::cout << "Initialising C++ to Mathematica link with model: " + frmodel << std::endl;
    
    // Declare FeynRules model class
    FeynRules model;
    
    // Open link to Mathematica
    model.create_wstp_link();
    
    // Load FeynRules
    model.load_feynrules();
    
    // Set the FeynRules model and load it up
    model.set_name(frmodel);
    model.load_model(frmodel);
    
    // Diagnositics -- check it is hermitian
    model.check_herm();
    
    
    // All done. Close the Mathematica link.
    model.close_wstp_link();
}
