#include <iostream>

#include "options.hpp"

void usage(std::string name)
{
    std::cerr << "Usage: " << name << " <option(s)>:\n\n"
              << "\t-h,--help\t\tShow this help message.\n"
              << "\t-p,--package PACKAGE\tSpecify the Mathematica package to use, 'feynrules' or 'sarah'.\n"
              << "\t-m,--model MODEL\tSpecify the name of the model in a given package.\n"
              << "\t-r,--restriction RESTRICTION\tSpecify the name of a restriction file for a given package.\n"
              << "\n"
              << "The user must specify at the very minimum -m and -p.\n"
              << "e.g. " << name << " -p feynrules -m SingletDM\n"
              << std::endl;
}

Options parse(int argc, char** argv, Options options)
{
    if (argc < 2) 
    {
        throw("");
    }
    for (int i=1; i<argc; ++i)
    {
        std::string arg = argv[i];
        if ((arg == "-h") || (arg == "--help")) 
        {
            throw("");
        } 
        else if ((arg == "-p") || (arg == "--package")) 
        {

            if (i+1 < argc) 
            { 
                options.package = argv[++i];
            } 
            else 
            { 
                throw("ERROR: -p option requires an argument.");
            }  
        } 
        else if ((arg == "-m") || (arg == "--model")) 
        {
            if (i+1 < argc) 
            { 
                options.model = argv[++i];
            } 
            else 
            { 
                throw("ERROR: -m option requires an argument.");
            }  
        }
        else if ((arg == "-r") || (arg == "--restriction")) 
        {
            if (i+1 < argc) 
            { 
                options.restriction = argv[++i];
            } 
            else 
            { 
                throw("ERROR: -r option requires an argument.");
            }  
        }  
    }
    if ((options.package != "sarah") && (options.package != "feynrules"))
    {
        throw("ERROR: no Mathematica package specified.");
    }
    if ((options.model.empty()))
    {
        throw("ERROR: no model file specified.");
    }
    return options;
}
