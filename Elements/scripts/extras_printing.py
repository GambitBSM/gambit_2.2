#!/usr/bin/env python
#
# GAMBIT: Global and Modular BSM Inference Tool
#*********************************************
# \file
#
#  Script to edit files, usually to make
#  sure no printing is used for standalones.
#  Also runs in reverse when not running standalone
#  to rewrite the printing. Extra files can be edited
#  provided they are also added as a dependency in the
#  utilities.cmake function extras_printing()
#
#*********************************************
#
#  Authors (add name and date if you modify):
#
#  \author Christopher Chang
#          (christopher.chang@uqconnect.edu.au)
#    \date 2021 Feb
#
#*********************************************

import pickle,sys

def Suspicious_Points(isStandalone):
    # Ensure suspicious_points.cpp doesn't print in standalone mode.

    includes = "\
#include \"gambit/Printers/printermanager.hpp\"\n\
#include \"gambit/Printers/baseprinter.hpp\"\n"

    printing = "\n\
    void Suspicious_point_exception::raise(const std::string& msg,int mycode, bool debug)\n\
    {\n\
      // get the printer pointer\n\
      Printers::BaseBasePrinter& printer = *(get_global_printer_manager()->printerptr);\n\
\n\
      printer.print(mycode, \"Suspicious Point Code\", Printers::get_main_param_id(\"Suspicious Point Code\"), printer.getRank(), Printers::get_point_id());\n"

    function = "\n\
    void Suspicious_point_exception::raise(const std::string& msg,int, bool debug)\n\
    {"

    filename = "./Elements/src/suspicious_points.cpp"

    with open(filename,"rt") as f:
        contents = f.read()

    if isStandalone:
        contents = contents.replace(includes,'// printer includes go here')
        contents = contents.replace(printing,function)

    else:
        contents = contents.replace('// printer includes go here',includes)
        contents = contents.replace(function,printing)

    with open(filename,"wt") as f:
        f.write(contents)
    return()


def main(isStandalone):

    # Add in further functions for other required edits
    Suspicious_Points(isStandalone)

if __name__ == "__main__":
   main(int(sys.argv[1]))
