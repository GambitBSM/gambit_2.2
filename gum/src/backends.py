"""
Master file containing all routines for modifying Backend interfaces.
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

def write_backend_patch(model, pristine_dir, patched_dir, backend, version):
    import subprocess
    outdir = "Outputs/"+model+"/Backends/patches/"+backend+"/"+version
    mkdir_if_absent(outdir)
    outfile = outdir+"/patch_"+backend+"_"+version+".dif"
    subprocess.call("diff -rupN "+pristine_dir+" "+patched_dir+" > "+outfile, shell=True)

def fix_pythia_lib(model, patched_dir):

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
                for x in processes:
                    f_new.write("  if ("+model+" || settings.flag(\""+model+":"+x[0]+"2"+x[1]+"\")) {\n")
                    f_new.write("    sigmaPtr = new Sigma_"+model+"_"+x[0]+"_"+x[1]+"();\n")
                    f_new.write("    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );\n")
                    f_new.write("  }\n")
    os.remove(old)
    os.rename(tmp, old)
