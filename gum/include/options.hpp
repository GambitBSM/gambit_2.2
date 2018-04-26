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

#endif
