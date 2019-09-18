#!/usr/bin/env python
#
#  GUM: GAMBIT Universal Models
#  ****************************
#  \file
#
#  Master file containing all routines for modifying Backend interfaces.
#
#  *************************************
#
#  \author Sanjay Bloor
#          (sanjay.bloor12@imperial.ac.uk)
#  \date 2017, 2018, 2019
#
#  \author Tomas Gonzalo
#          (tomas.gonzalo@monash.edu)
#  \date 2019 July, August
#
#  **************************************


import os
import re

from setup import *
from files import *
from parse import *
from cmake_variables import *

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

    if outputs.vev:

        if os.path.exists("./../Backends/installed/vevacious/VevaciousPlusPlus/1.0/ModelFiles/"):
            print("Found vevacious.")
        else:
            raise GumError(("\n\nNo vevacious installation found. Please go to into"
                            " the GAMBIT build directory and do"
                            ":\n   make vevacious"))


    print("\nAll backends found -- connecting to Mathematica!\n")


def write_backend_patch(output_dir, pristine_dir, patched_dir, backend, version):
    """
    Writes a backend patch (i.e. just a diff) from a pristine 
    and a modified version of a backend.
    """
    import subprocess
    full_output_dir = output_dir+"/Backends/patches/"+backend+"/"+version
    mkdir_if_absent(full_output_dir)
    outfile = full_output_dir+"/patch_"+backend+"_"+version+".dif"
    pristine_parts = os.path.split(pristine_dir)
    cwd = os.getcwd()
    os.chdir(pristine_parts[0])
    subprocess.call("diff -rupN "+pristine_parts[1]+" "+patched_dir+" > "+outfile, shell=True)
    os.chdir(cwd)


def write_new_default_bossed_version(backend, version, output_dir):

    import re

    # The path to the original file in GAMBIT
    path = "/Backends/include/gambit/Backends/"
    filename = "default_bossed_versions.hpp"
    old = ".."+path+filename

    # Sort out the path to the candidate replacement
    newdir = output_dir+path
    mkdir_if_absent(newdir)
    new = newdir+filename

    # The signature of the line we want to add/replace
    signature = "#define  Default_"+backend+" "

    # Flag indicating that GUM section exists already
    comment_exists = False

    # Work through the old version of the file and add/replace this entry
    with open(old) as f_old, open(new, 'w') as f_new:
        for line in f_old:
            if not signature in line: f_new.write(line)
            if "// Defaults added by GUM" in line:
                if not line.endswith("\n"): f_new.write("\n")
                comment_exists = True
        if not comment_exists:
            f_new.write("\n// Defaults added by GUM (do not remove this comment).\n")
        f_new.write(signature+re.sub("\.", "_", version)+"\n")


def add_to_backend_locations(backend_name, backend_location, version_number, reset_dict):
    """
    Adds an entry to backend_locations.yaml for a new backend.
    """

    # Check to see if backend_locations.yaml exists; if not then we'll use
    # backend_locations.yaml.default (as long as it hasn't been removed)

    if not os.path.isfile("./../config/backend_locations.yaml.default"):
        raise GumError("backend_locations.yaml.default is missing. What have you done to GAMBIT!?")

    target = "backend_locations.yaml"

    if not os.path.isfile("./../config/"+target):
        target = "backend_locations.yaml.default"

    # Add the new backend before the examples stuff.
    linenum = 0
    with open("./../config/"+target) as f:
      for num, line in enumerate(f, 1):
          if "Example" in line:
              linenum = num
              break

    contents = ("# Added by GUM\n"
                "{0}:\n"
                "  {1}:         ./Backends/installed/{2}"
                "\n"
                "\n"
                ).format(backend_name, version_number, backend_location)

    # Write the changes
    amend_file(target, "config", contents, linenum-1, reset_dict)


def add_to_backends_cmake(contents, reset_dict, linenum=0, string_to_find=""):
    """
    Adds an entry to backends.cmake, either with a line number
    or with a string to find, then write afterwards.
    """

    if ((linenum == 0) and (string_to_find == "")):
        raise GumError(("\n\tYou need to pass either a line number, or a string for "
                        "me to find, if you want to amend backends.cmake."))

    # If the user specifies a 
    if (linenum != 0):
        amend_file("backends.cmake", "cmake", contents, linenum-1, reset_dict)

    else:
        present, linenum = find_string("backends.cmake", "cmake", string)
        if present: amend_file("backends.cmake", "cmake", contents, linenum-1, reset_dict)


def write_backend_frontends(output_dir, model_name, backend_name, backend_version, header_content, src_content) :
    """
    Adds frontend headers and source files for backend to output directory
    """

    # Header and source directories
    header_dir = output_dir + "/Backends/include/gambit/Backends/frontends/"
    src_dir = output_dir + "/Backends/src/frontends/"

    mkdir_if_absent(header_dir)
    mkdir_if_absent(src_dir)

    # Name of files (except extension)
    clean_model_name = model_name.replace("-","")
    safe_version = backend_version.replace(".","_")
    filename = backend_name + "_" + clean_model_name + "_" + safe_version

    # Write to files now
    with open(header_dir + "/" + filename + ".hpp", 'w') as f_header :
      f_header.write(header_content)

    with open(src_dir + "/" + filename + ".cpp", 'w') as f_src :
      f_src.write(src_content)
