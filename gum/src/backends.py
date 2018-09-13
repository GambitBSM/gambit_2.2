"""
Master file containing all routines for modifying Backend interfaces.
Think CalcHEP, MadGraph, SPheno.
"""

import os

from setup import *
from files import *
from parse import *

def check_backends(outputs):
    """
    Diagonostics to check all backends exist in the GAMBIT repository.
    """
    ## TO DO - SPheno, MadGraph, FlexibleSUSY, Vevacious...

    if not isinstance(outputs, Outputs):
        raise GumError("\nRequested output not passed as class Outputs.\n")

    print("\nChecking for backends before we get going...\n")

    # CalcHEP
    if outputs.ch:

      if os.path.exists("./../Backends/installed/calchep/3.6.27/models/"):
          print("Found CalcHEP.")
      else:
          raise GumError(("\n\nNo CalcHEP installation found. Please go to into"
                          " the GAMBIT build directory and do"
                          ":\n   make calchep"))

    print("\nAll backends found -- connecting to Mathematica!\n")

def add_calchep_switch(model_name, spectrum):
    """
    Adds an 'if ModelInUse()' switch to the CalcHEP frontend to make GAMBIT
    point to the correct CalcHEP files.
    """

    # Scan-level
    src_sl = dumb_indent(2, (
        "if (ModelInUse(\"{0}\"))\n"
        "{{\n"
        "BEpath = backendDir + \"/../models/{0};\n"
        "path = BEpath.c_str();\n"
        "modeltoset = (char*)malloc(strlen(path)+11);\n"
        "sprintf(modeltoset, \"%s\", path);\n"
        "}}\n\n"
    ).format(model_name))

    # Point-level
    src_pl = dumb_indent(2, (
           "if (ModelInUse(\"{0}\"))\n"
           "{{\n"
           "// Obtain model contents\n"
           "static const SpectrumContents::{0} {0}_contents;\n\n"
           "// Obtain list of all parameters within model\n"
           "static const std::vector<SpectrumParameter> {0}_params = "
           "{0}_contents.all_parameters();\n\n"
           "// Obtain spectrum information to pass to CalcHEP\n"
           "const Spectrum& spec = *Dep::{1};\n\n"
           "Assign_All_Values(spec, {0}_params);\n"
           "}}\n\n"
    ).format(model_name, spectrum))

    # to do -- also ADD_MODEL()
    header = (
           "BE_INI_CONDITIONAL_DEPENDENCY({0}, Spectrum, {1})\n"
    ).format(spectrum, model_name)

    return indent(src_sl), indent(src_pl), header




