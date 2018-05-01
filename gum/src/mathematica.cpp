/// Contains all routines for interface between C++ and Mathematica packages.

#include <iostream>

#include "options.hpp"
#include "feynrules.hpp"
#include "sarah.hpp"

int main(int argc, char** argv)
{
    
    // Initialise options for input.
    std::vector<Particle> partlist;
    std::vector<Parameter> paramlist;

    // Attempt to parse the user's command line input...
    try
    {
        Options options = parse(argc, argv);
            
        // Pick FeynRules or SARAH, then crack on.
        if (options.package() == "feynrules")
        {
            #ifndef HAVE_FEYNRULES
            std::cerr << "\n ERROR: You can't call FeynRules if you haven't linked to it!!" 
                      << "\n  Try 'cmake -D FR=path_to_feynrules ..' in your build directory.\n\n";
            return 0;
            #endif
            
            all_feynrules(options, partlist, paramlist);
        }
        else if (options.package() == "sarah")
        {
        
            #ifndef HAVE_SARAH
            std::cerr << "\n ERROR: You can't call SARAH if you haven't linked to it!!" 
                      << "\n  Try 'cmake -D SARAH=path_to_sarah ..' in your build directory.\n\n";
            return 0;
            #endif
        
            all_sarah(options);
        }
    }
    catch(const char* e)
    {
        std::cerr << e << std::endl << std::endl;
        usage(argv[0]);
        return 0;
    }

    std::cout << "C++ to Mathematica interface finished." << std::endl;
        
    return 0;
}
