#ifndef OPTIONS
#define OPTIONS

#include <iostream>

class Options
{

    std::string optpackage;
    std::string optmodel;
    std::string optrestriction;
    std::string optLTot;

    public:
        Options(std::string package, std::string model, std::string restriction, std::string lagrangian = "LTotal")
        {
            optpackage = package;
            optmodel = model;
            optrestriction = restriction;
            optLTot = lagrangian;
        }

        std::string model() { return optmodel; }
        std::string package() { return optpackage; }
        std::string restriction() { return optrestriction; }
        std::string lagrangian() { return optLTot; }

};


class Outputs
{
    std::string ch;
    std::string mg;

    public:

        std::string get_ch() { return ch; }
        std::string get_mg() { return mg; }

        void set_ch(std::string chdir) { ch = chdir; }
        void set_mg(std::string mgdir) { mg = mgdir; }
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
    bool selfconj;
    std::string antipartname;

    public:

        // Needed for Boost.python interface
        bool operator==(const Particle&) {return false;}
        bool operator!=(const Particle&) {return true;}

        Particle(int pdg, std::string name, int spinX2, std::string full_name, bool SM, std::string mass, std::string antiname)
        {
            pdgcode = pdg;
            partname = name;
            doublespin = spinX2;
            fullname = full_name;
            standardmodel = SM;
            partmass = mass;
            antipartname = antiname;

            if (name == antiname)
            {
                selfconj = true;
            }
            else
            {
                selfconj = false;
            }

        }

        int pdg() { return pdgcode; }
        std::string name() { return partname; }
        bool SM() { return standardmodel; }
        int spinX2() { return doublespin; }
        std::string mass() { return partmass; }
        bool SC() { return selfconj; }
        std::string antiname() { return antipartname; }


        // Unsure if needed yet.
        int colourX3;

};

class Parameter
{

    std::string paramname;
    std::string blockname;

    public:

        // Needed for Boost.python interface
        bool operator==(const Parameter&) {return false;}
        bool operator!=(const Parameter&) {return true;}

        Parameter(std::string name, std::string block)
        {
            paramname = name;
            blockname = block;
        }

        std::string name() { return paramname; }
        std::string block() { return blockname; }

};

#endif
