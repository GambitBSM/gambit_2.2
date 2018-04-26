/// Contains all routines for interface between C++ and Mathematica packages.

#include <iostream>

#include "options.hpp"
#include "feynrules.hpp"
#include "sarah.hpp"

int main(int argc, char** argv)
{
    
    // Initialise options for input.
    Options options;

    // Attempt to parse the user's command line input...
    try
    {
        options = parse(argc, argv, options);
    }
    catch(const char* e)
    {
        std::cerr << e << std::endl << std::endl;
        usage(argv[0]);
        return 0;
    }
    
    // Pick FeynRules or SARAH, then crack on.
    if (options.package == "feynrules")
    {
        all_feynrules(options);
    }
    else if (options.package == "sarah")
    {
        all_sarah(options);
    }
    
    std::cout << "C++ to Mathematica interface: done." << std::endl;
        
    return 0;
}
