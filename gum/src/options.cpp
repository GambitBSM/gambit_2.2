#include <iostream>

#include "options.hpp"

void usage(std::string name)
{
    std::cerr << "Usage: " << name << " <option(s)>:\n\n"
              << "\t-h,--help\t\t\tShow this helpful help message.\n"
              << "\t-p,--package PACKAGE\t\tSpecify the Mathematica package to use, 'feynrules' or 'sarah'.\n"
              << "\t-m,--model MODEL\t\tSpecify the name of the model in a given package. This should live in $GUM/Models/<modelname>/<modelname>.[m/fr]\n"
              << "\n\t ** FeynRules specific ** \n\n"
              << "\t-b,--basemodel BASEMODEL\t\tSpecify the base model to use, if the FeynRules model file depends on another (e.g. the SM). This should live in $GUM/Models/<modelname>/<modelname>.[m/fr]\n"
              << "\t-r,--restriction RESTRICTION\tSpecify the name of a restriction file for FeynRules.\n"
              << "\t-L,--lagrangian LAGRANGIAN\tSpecify the total Lagrangian for a FeynRules model.\n"
              << "\n"
              << "The user must specify at the very minimum -m and -p.\n"
              << "e.g. " << name << " -p feynrules -m SingletDM\n"
              << std::endl;
}

Options parse(int argc, char** argv)
{
    std::string package;
    std::string model;
    std::string basemodel;
    std::string restriction;
    std::string lagrangian;

    if (argc < 2)
    {
        throw("Here are some tips to get you going.");
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
                package = argv[++i];
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
                model = argv[++i];
            }
            else
            {
                throw("ERROR: -m option requires an argument.");
            }
        }        
        else if ((arg == "-b") || (arg == "--basemodel"))
        {
            if (i+1 < argc)
            {
                basemodel = argv[++i];
            }
            else
            {
                throw("ERROR: -b option requires an argument.");
            }
        }
        else if ((arg == "-r") || (arg == "--restriction"))
        {
            if (i+1 < argc)
            {
                restriction = argv[++i];
            }
            else
            {
                throw("ERROR: -r option requires an argument.");
            }
        }
        else if ((arg == "-L") || (arg == "--Lagrangian"))
        {
            if (i+1 < argc)
            {
                lagrangian = argv[++i];
            }
            else
            {
                throw("ERROR: -L option requires an argument.");
            }
        }
    }
    if ((package != "sarah") && (package != "feynrules"))
    {
        throw("ERROR: no Mathematica package specified.");
    }
    if ((model.empty()))
    {
        throw("ERROR: no model file specified.");
    }
    if ((package == "sarah") && (not restriction.empty()))
    {
        throw("ERROR: restriction file added for SARAH, only works for FeynRules.");
    }
    if ((package == "sarah") && (not basemodel.empty()))
    {
        throw("ERROR: base model specified for SARAH, only works for FeynRules.");
    }
    Options options(package, model, basemodel, restriction, lagrangian);
    return options;
}
