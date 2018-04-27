#ifndef OPTIONS
#define OPTIONS

#include <iostream>

class Options
{
    public:
        std::string package;
        std::string model;
        std::string restriction;

};

void usage(std::string);

Options parse(int, char**, Options);


class Particle
{

    int pdgcode;
    std::string partname;
    int doublespin;
    std::string fullname;
    bool standardmodel;
    
    public:
        Particle(int pdg, std::string name, int spinX2, std::string full_name, bool SM)
        {
            pdgcode = pdg;
            partname = name;
            doublespin = spinX2;
            fullname = full_name;
            standardmodel = SM;
        }
        
        int pdg() { return pdgcode; }
        std::string name() { return partname; }
        bool SM() { return standardmodel; }
        int spinX2() { return doublespin; }

        
        // Unsure if needed yet. 
        int colourX3;
        
};

class Parameter
{

    std::string paramname;
    bool standardmodel;
        
    public:
        Parameter(std::string name, bool SM)
        {
            paramname = name;
            standardmodel = SM;
        }
        
        std::string name() { return paramname; }
        bool SM() { return standardmodel; }
        
};

#endif
