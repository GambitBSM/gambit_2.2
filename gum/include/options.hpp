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
    std::string partmass;
    
    public:
        Particle(int pdg, std::string name, int spinX2, std::string full_name, bool SM, std::string mass)
        {
            pdgcode = pdg;
            partname = name;
            doublespin = spinX2;
            fullname = full_name;
            standardmodel = SM;
            partmass = mass;
        }
        
        int pdg() { return pdgcode; }
        std::string name() { return partname; }
        bool SM() { return standardmodel; }
        int spinX2() { return doublespin; }
        std::string mass() { return partmass; }

        
        // Unsure if needed yet. 
        int colourX3;
        
};

class Parameter
{

    std::string paramname;
    std::string blockname;
        
    public:
        Parameter(std::string name, std::string block)
        {
            paramname = name;
            blockname = block;
        }
        
        std::string name() { return paramname; }
        std::string block() { return blockname; }
        
};

#endif
