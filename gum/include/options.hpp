#ifndef OPTIONS
#define OPTIONS

#include <iostream>

class Options
{

    std::string optpackage;
    std::string optmodel;
    std::string optrestriction;
    
    public:
        Options(std::string package, std::string model, std::string restriction)
        {
            optpackage = package;
            optmodel = model;
            optrestriction = restriction;
        }
        
        std::string model() { return optmodel; }
        std::string package() { return optpackage; }
        std::string restriction() { return optrestriction; }

};

void usage(std::string);

Options parse(int, char**);


class Particle
{

    int pdgcode;
    std::string partname;
    int doublespin;
    std::string fullname;
    bool standardmodel;
    std::string partmass;
    
    public:
        
        bool operator==(const Particle& other) {return false;}
        bool operator!=(const Particle& other) {return true;}
    
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
    
        bool operator==(const Parameter& other) {return false;}
        bool operator!=(const Parameter& other) {return true;}
        
        Parameter(std::string name, std::string block)
        {
            paramname = name;
            blockname = block;
        }
        
        std::string name() { return paramname; }
        std::string block() { return blockname; }
        
};

#endif
