"""
Master file containing all routines for modifying Backend interfaces.
"""

import os

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

    print("\nAll backends found -- connecting to Mathematica!\n")

def add_calchep_switch(model_name, spectrum):
    """
    Adds an 'if ModelInUse()' switch to the CalcHEP frontend to make GAMBIT
    point to the correct CalcHEP files.
    """

    # Scan-level
    src_sl = dumb_indent(4, (
        "if (ModelInUse(\"{0}\"))\n"
        "{{\n"
        "BEpath = backendDir + \"/../models/{0}\";\n"
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

    # to do -- also ALLOW_MODEL()
    header = (
           "BE_INI_CONDITIONAL_DEPENDENCY({0}, Spectrum, {1})\n"
    ).format(spectrum, model_name)

    return indent(src_sl), indent(src_pl), header

def write_backend_patch(output_dir, pristine_dir, patched_dir, backend, version):
    import subprocess
    full_output_dir = output_dir+"/Backends/patches/"+backend+"/"+version
    mkdir_if_absent(full_output_dir)
    outfile = full_output_dir+"/patch_"+backend+"_"+version+".dif"
    pristine_parts = os.path.split(pristine_dir)
    cwd = os.getcwd()
    os.chdir(pristine_parts[0])
    subprocess.call("diff -rupN "+pristine_parts[1]+" "+patched_dir+" > "+outfile, shell=True)
    os.chdir(cwd)

def fix_pythia_lib(model, patched_dir, pythia_groups):
    """
    Routine to patch the new Pythia - adding new matrix elements to the
    Process Container -- and the shared library too.
    """

    import re

    # Move the matrix element sources and headers to where they will get compiled into the Pythia library
    process_dir = patched_dir+"/Processes_"+model+"/"
    source_dir = patched_dir+"/src/"
    inc_dir = patched_dir+"/include/"
    files = os.listdir(process_dir)
    headers = [x for x in files if x.endswith(".h")]
    sources = [x for x in files if x.endswith(".cc")]
    for x in headers:
        os.rename(process_dir+"/"+x, inc_dir+"/"+x)
    for x in sources:
        os.rename(process_dir+"/"+x, source_dir+"/"+x)

    # Scrape the list of processes from the header names
    processes = []
    for x in headers:
        if (x.startswith("Sigma_")):
            states = re.split("_", re.sub("Sigma_"+model+"_", "", re.sub("\.h", "", x)))
            with open(inc_dir+"/"+x) as f:
                for line in f:
                    if "// Process: " in line:
                        process_string = re.sub(" WEIGHTED.*\n", "", re.sub("// Process: ", "", line))
                        process_string = re.sub(">", "&rarr;", process_string)
                        states.append("<ei>"+process_string+"<\ei>")
            processes.append(states)

    # Get rid of the leftover process directory made by MadGraph
    remove_tree_quietly(process_dir)

    # Write the xml doc file for the new processes
    with open(patched_dir+"/share/Pythia8/xmldoc/"+model+"Processes.xml", 'w') as f:
        f.write("<flag name=\""+model+":all\" default=\"off\">\n")
        f.write("Common switch for production of "+model+" processes. Added by GAMBIT.\n")
        f.write("</flag>\n")

        # Go through pythia_groups to add each individual flag, to
        # select groups of subprocesses
        if pythia_groups:
            for group in pythia_groups:
                for k, v in group.items():
                    f.write("<flag name=\""+model+k+":all\" default=\"off\">\n")
                    f.write("Common switch for production of "+model+" processes, involving the group of particles ["+', '.join(v)+"] as external legs ONLY. Added by GAMBIT.\n")
                    f.write("</flag>\n")

        # Invidiual processes
        for x in processes:
            f.write("\n")
            f.write("<flag name=\""+model+":"+x[0]+"2"+x[1]+"\" default=\"off\">\n")
            f.write("Switch for "+model+" process "+x[2]+". Added by GAMBIT.\n")
            f.write("</flag>\n")

    # Add the new processes to the master Pythia index.xml file
    old = patched_dir+"/share/Pythia8/xmldoc/Index.xml"
    tmp = old+".temp"
    with open(old) as f_old, open(tmp, 'w') as f_new:
        for line in f_old:
            f_new.write(line)
            if "<aidx href=\"SUSYProcesses\">SUSY</aidx><br/>" in line:
                f_new.write("&nbsp;&nbsp;--&nbsp;&nbsp;\n")
                f_new.write("<aidx href=\""+model+"Processes\">"+model+"</aidx><br/>\n")
    os.remove(old)
    os.rename(tmp, old)

    # Add the new processes to the Pythia Process Container
    old = source_dir+"ProcessContainer.cc"
    tmp = old+".temp"
    with open(old) as f_old, open(tmp, 'w') as f_new:
        for line in f_old:
            f_new.write(line)
            if "#include \"Pythia8/SigmaSUSY.h\"" in line:
                for x in headers:
                    if x.startswith("Sigma"):
                        f_new.write("#include \""+x+"\"\n")

            if "// End of SUSY processes." in line:

                f_new.write("\n")
                f_new.write("  // Set up requested objects for "+model+" processes. Added by GAMBIT.\n")
                f_new.write("  bool "+model+" = settings.flag(\""+model+":all\");\n")

                # Go through pythia_groups to add each individual flag, to
                # select groups of subprocesses
                if pythia_groups:
                    for group in pythia_groups:
                        for k, v in group.items():
                            f_new.write("  bool "+model+k+" = settings.flag(\""+model+k+":all\");\n")

                # Add each process
                for x in processes:

                    process_syntax = x[2]
                    # Get the in and out states
                    in_states =  re.search('<ei>(.*)&rarr', process_syntax).group(1).split()
                    out_states = re.search(';(.*)<\\\\ei>', process_syntax).group(1).split()
                    external_states = in_states+out_states

                    # If any of the particles belong to any Pythia group,
                    # add them as a conditional to be selectable from
                    # the yaml file in a GAMBIT scan.

                    switches = ""
                    if pythia_groups:
                        for group in pythia_groups:
                            for k, v in group.items():
                                vtemp = [i.lower() for i in v] # Pythia makes everything lowercase
                                if any([i.lower() in vtemp for i in external_states]):
                                    switches += " {0} ||".format(model+k)

                    f_new.write("  if ("+model+" ||{0} settings.flag(\"".format(switches)+model+":"+x[0]+"2"+x[1]+"\")) {\n")
                    f_new.write("    sigmaPtr = new Sigma_"+model+"_"+x[0]+"_"+x[1]+"();\n")
                    f_new.write("    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );\n")
                    f_new.write("  }\n")
    os.remove(old)
    os.rename(tmp, old)

def write_boss_config_for_pythia(model, output_dir):

    # Sort out the paths
    path = "/Backends/scripts/BOSS/configs"
    filename = "/pythia_"+model.lower()+"_8_"+base_pythia_version+".py"
    full_output_dir = output_dir+path
    mkdir_if_absent(full_output_dir)
    template = ".."+path+"/pythia_8_"+base_pythia_version+".py"
    outfile = full_output_dir+filename

    # Write the actual BOSS config file
    with open(outfile, 'w') as f_new, open(template) as f_old:
        for line in f_old:
          if "gambit_backend_name    = 'Pythia'" in line:
              f_new.write("gambit_backend_name    = 'Pythia_"+model+"'\n")
          else:
              f_new.write(line)
          if "#  Configuration module for BOSS  #" in line:
              f_new.write("#  ----brought to you by GUM----  #\n")


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


def add_new_pythia_to_backends_cmake(model, output_dir):

    # The string that will commence the block to be added by GUM
    signature = "# Pythia with matrix elements for "+model+" (brought to you today by the letters G, U and M)."

    # The path to the original file in GAMBIT
    old = "../cmake/backends.cmake"
    # Sort out the path to the candidate replacement
    newdir = output_dir+"/cmake"
    mkdir_if_absent(newdir)
    new = newdir+"/backends.cmake"

    # Initialise flags to indicate place in the original file
    in_duplicate = False
    passed_pythia = False
    wrote_entry = False

    # Open old and new files and iterate through the old one, writing to the new as we go.
    with open(old) as f_old, open(new, 'w') as f_new:
        for line in f_old:
            # Have we spotted a previous modification by GUM?  If so, overwrite it.
            if signature in line: in_duplicate = True
            # Have we spotted the vanilla pythia entry yet?
            if "set(name \"pythia\")" in line: passed_pythia = True
            if not in_duplicate: f_new.write(line)
            if not wrote_entry and passed_pythia and "set_as_default_version(\"backend\" ${name} ${ver})" in line:
                to_write = "endif()\n"\
                           "\n"+signature+"\n"\
                           "set(model \""+model.lower()+"\")\n"\
                           "set(name \"pythia_${model}\")\n"\
                           "set(ver \"8."+base_pythia_version+"\")\n"\
                           "set(lib \"libpythia8\")\n"\
                           "set(dl \"http://home.thep.lu.se/~torbjorn/pythia8/pythia8"+base_pythia_version+".tgz\")\n"\
                           "set(md5 \""+pythia_md5+"\")\n"\
                           "set(dir \"${PROJECT_SOURCE_DIR}/Backends/installed/${name}/${ver}\")\n"\
                           "set(patch1 \"${PROJECT_SOURCE_DIR}/Backends/patches/${name}/${ver}/patch_${name}_${ver}.dif\")\n"\
                           "set(patch2 \"${PROJECT_SOURCE_DIR}/Backends/patches/pythia/${ver}/patch_pythia_${ver}.dif\")\n"\
                           "check_ditch_status(${name} ${ver})\n"\
                           "if(NOT ditched_${name}_${ver})\n"\
                           "  ExternalProject_Add(${name}_${ver}\n"\
                           "    DOWNLOAD_COMMAND ${DL_BACKEND} ${dl} ${md5} ${dir} ${name} ${ver}\n"\
                           "    SOURCE_DIR ${dir}\n"\
                           "    BUILD_IN_SOURCE 1\n"\
                           "    PATCH_COMMAND patch -p1 < ${patch1}\n"\
                           "          COMMAND patch -p1 < ${patch2}\n"\
                           "    CONFIGURE_COMMAND ./configure --enable-shared --cxx=\"${CMAKE_CXX_COMPILER}\" --cxx-common=\"${pythia_CXXFLAGS}\" --cxx-shared=\"${pythia_CXX_SHARED_FLAGS}\" --lib-suffix=\".so\"\n"\
                           "    BUILD_COMMAND ${CMAKE_MAKE_PROGRAM} CXX=\"${CMAKE_CXX_COMPILER}\" lib/${lib}.so\n"\
                           "    INSTALL_COMMAND \"\"\n"\
                           "  )\n"\
                           "  BOSS_backend(${name} ${ver})\n"\
                           "  add_extra_targets(\"backend\" ${name} ${ver} ${dir} ${dl} distclean)\n"\
                           "  set_as_default_version(\"backend\" ${name} ${ver})\n"
                f_new.write(to_write)
                wrote_entry = True
            # We've reached the end of the previous modification by GUM, so remove the hold on repeating lines from the old file.
            if in_duplicate and "endif()" in line: in_duplicate = False


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

    contents = ("#Added by GUM\n"
                "{0}:\n"
                "  {1}:         ./Backends/installed/{2}"
                "\n"
                "\n"
                ).format(backend_name, version_number, backend_location)

    # Write the changes
    amend_file(target, "config", contents, linenum-1, reset_dict)


