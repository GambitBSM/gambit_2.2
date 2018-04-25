"""
Master module for Mathematica related routines.
"""

import os

from setup import *

def find_mathematica():
    """
    Uses gambit cmake variables to ensure C++ to Mathematica pathway
    is established.
    """
    
    print("Trying to find Mathematica...")
    
    cmake_vars = "../cmake/include/gambit/cmake/cmake_variables.hpp"
    lookup = "#define HAVE_MATHEMATICA 1"
    
    if not os.path.exists(cmake_vars):
        msg = (
            "\n"
            "\n"
            "File {0} not found. Please do the following:\n\n"
            "  cd ..\n"
            "  mkdir build\n"
            "  cd build\n"
            "  cmake ..\n"
        ).format(cmake_vars)
        raise GumError(msg)
    
    if lookup in open(cmake_vars).read():
        print("Mathematica found.")
        return True
        
    return False    
    
    


    
def find_sarah():
    """
    Tries to find SARAH.
    """

def find_feynrules():
    """
    Tries to find FeynRules.
    """
