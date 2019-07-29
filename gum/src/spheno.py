#!/usr/bin/env python
#
#  GUM: GAMBIT Universal Models
#  ****************************
#  \file
#
#  Routines for SPheno output from SARAH
#
#  *************************************
#
#  \author Sanjay Bloor
#          (sanjay.bloor12@imperial.ac.uk)
#  \date 2019 July
#
#  \author Tomas Gonzalo
#          (tomas.gonzalo@monash.edu)
#  \date 2019 July
#
#  **************************************

from files import *
from cmake_variables import *
from distutils.dir_util import copy_tree
from collections import defaultdict
import re
import os.path

# An empty defaultdict in scope, so we don't save to the 
# mug file, but can still use the writing routines in 
# the 'files' module.
d = defaultdict(list)

def copy_spheno_files(model_name, output_dir, spheno_oob_path, sarah_spheno_path):
    """
    Creates a copy of SPheno output in the 
    Outputs/... folder.
    Then create another one, for the patched version.
    """

    dirs = [output_dir + "/SPheno", output_dir+"/SPheno_patched"]

    for newdir in dirs:

        # Remove the directory if it already exists, then make it
        remove_tree_quietly(newdir)
        mkdir_if_absent(newdir)

        # Copy SPheno to the Output directory
        copy_tree(spheno_oob_path, newdir)

        # Now copy from SARAH to the Output directory
        modeldir = newdir + "/" + model_name
        mkdir_if_absent(modeldir)
        copy_tree(sarah_spheno_path, modeldir)

    print("SPheno files moved to output directory, and a copy made.")



"""
PATCHING
"""

def patch_spheno(model_name, patch_dir):
    """
    Applies all patches to SPheno in the GUM 
    Outputs/... directory.
    """

    patch_spheno_makefile(model_name, patch_dir)
    patch_spheno_model_makefile(model_name, patch_dir)
    patch_spheno_src_makefile(model_name, patch_dir)
    patch_control(model_name, patch_dir)
    patch_brs(model_name, patch_dir)
    patch_loopfunctions(model_name, patch_dir)
    patch_spheno_model(model_name, patch_dir)

    if model_name == "MSSM" or model_name == "NMSSM" :
        # TODO: if gum.is_susy: ...
        patch_model_data(model_name, patch_dir)
        patch_3_body_decays_susy(model_name, patch_dir)

    print("SPheno files patched.")


def patch_spheno_makefile(model_name, patch_dir):
    """
    Patches $SPheno/Makefile
    """

    filename = patch_dir + "/Makefile"

    if not os.path.exists(filename):
        raise GumError(("Tried to find the file located at " + filename +
                        " but it does not seem to exist!"))

    with open(filename, 'w') as f :
        content = "# please put here your preferred F95/F2003 compiler\n"\
                  "# the options in src/Makefile have been put for the\n"\
                  "# cases NAG's nagfor, gfortran, g95, Lahey's lf95 and Intels ifort\n"\
                  "# Please uncomment the corresponding line\n"\
                  "# F90 = nagfor\n"\
                  "F90 = gfortran\n"\
                  "# F90 = g95\n"\
                  "# F90 = lf95\n"\
                  "# F90 = ifort\n"\
                  "Model = src\n"\
                  "version = 400.00\n"\
                  "all: bin/SPheno lib/libSPheno"+model_name+".so\n"\
                  "bin/SPheno:\n"\
                  "\tcd ${Model} ; ${MAKE} F90=${F90} version=${version}\n"\
                  "lib/libSPheno"+model_name+".so:\n"\
                  "\tcd ${Model} ; ${MAKE} $@ F90=${F90} version=${version}\n"\
                  "clean:\n"\
                  "\trm -f *.o *~ */*.o */*~\n"\
                  "cleanall:\n"\
                  "\trm -f bin/SPheno lib/*.a lib/*.so *.o *~ */*.o */*~ include/*\n"\
                  ".PHONY: bin/SPheno lib/libSPheno"+model_name+" clean cleanall"

        f.write(content)


def patch_spheno_model_makefile(model_name, patch_dir):
    """
    Patches $SPheno/<MODEL>/Makefile
    """
    
    filename = patch_dir + "/" + model_name + "/Makefile" 
    temp_filename = filename + "_temp"  

    if not os.path.exists(filename):
        raise GumError(("Tried to find the file located at " + filename +
                        " but it does not seem to exist!"))

    with open(filename, 'r') as f, open(temp_filename, 'w') as g :
        skip_next_lines = False
        for line in f :
            if line.startswith("name = ") :
                g.write(line)
                g.write("shared = lib/libSPheno" + model_name + ".so\n")
            elif line.startswith("F90=gfortran") :
                content = "ifeq (${F90},/usr/bin/gfortran)\n"\
                          "override F90=gfortran\n"\
                          "endif\n"\
                          "# Intels ifort,debug modus\n"\
                          "ifeq (${F90},ifort)\n"\
                          "F90=ifort\n"\
                          "comp= -c -g -fPIC -module ${Mdir} -I${InDir}\n"\
                          "LFlagsB= -g -fPIC\n"\
                          "endif\n"\
                          "# gfortran\n"\
                          "ifeq (${F90},gfortran)\n"\
                          "comp= -c -g -fPIC --free-line-length-none -J${Mdir} -I${InDir}\n"\
                          "LFlagsB= -g -fPIC\n"\
                          "endif\n"\
                          "# g95\n"\
                          "ifeq (${F90},g95)\n"\
                          "comp= -c -O -fPIC -fmod=${Mdir} -I${InDir}\n"\
                          "LFlagsB= -O -fPIC\n"\
                          "endif\n"\
                          "# Lahey F95 compiler\n"  \
                          "ifeq (${F90},lf95)\n"\
                          "comp=-c -O -fPIC -M ${Mdir} -I${InDir}\n"\
                          "LFlagsB=-O -fPIC\n"\
                          "endif\n"\
                          "# NAG f95/2003\n"\
                          "ifeq (${F90},nagfor)\n"\
                          "comp= -c -O -fPIC -mdir ${Mdir} -I${InDir}\n"\
                          "LFlagsB= -O -fPIC -DONLYDOUBLE -mdir ${MDir} -I${InDir}\n"\
                          "endif\n"\
                          ".SUFFIXES : .o .ps .f90 .F90 .a\n"
                g.write(content)
                skip_next_lines = True
            elif line.startswith("bin/SPheno" + model_name) :
                content = "${shared}:\n"\
                          "\tcd ../src ; ${MAKE} F90=${F90}\n"\
                          "\t${MAKE} F90=${F90} ${name}\n"\
                          "\t${MAKE} F90=${F90} SPheno" + model_name + ".o\n"\
                          "\t${F90} -c -fPIC TwoLoopMasses/effpotasat.f\n"\
                          "\t${F90} -shared -fPIC -o ../${shared} ${LFlagsB} SPheno" + model_name + ".o effpotasat.o ../lib/libSPheno" + model_name + ".a ../lib/libSPheno.a\n"\
                          "bin/SPheno" + model_name + ":\n"\
                          "ifeq (${cVersion},1)\n"\
                          "\tcd ../src ; ${MAKE} F90=${F90}\n"\
                          "\t${MAKE} F90=${F90} ${name}\n"\
                          "\t${MAKE} F90=${F90} SPheno" + model_name + ".o\n"\
                          "\t${F90} -c -fPIC TwoLoopMasses/effpotasat.f\n"\
                          "\t${F90} -o -fPIC SPheno" + model_name + " ${LFlagsB} SPheno" + model_name + ".o effpotasat.o ../lib/libSPheno" + model_name + ".a ../lib/libSPheno.a\n"\
                          "\tmv SPheno" + model_name + "../bin\n"\
                          "\trm SPheno" + model_name + ".o\n"
                g.write(content)
                for i in range(7): next(f)
                skip_next_lines = False
            elif line.startswith("cleanall:") :
                 g.write(line)
                 g.write("\trm -f bin/SPheno3 lib/*.a lib/*.so *~ */*.o */*~ include/*\n")
                 next(f)
            else :
                if not skip_next_lines :
                    g.write(line)

    os.remove(filename)
    os.rename(temp_filename, filename)

def patch_spheno_src_makefile(model_name, patch_dir):
    """
    Patches $SPheno/src/Makefile
    """

    filename = patch_dir + "/src/Makefile"
    temp_filename = filename + "_temp"

    if not os.path.exists(filename):
        raise GumError(("Tried to find the file located at " + filename +
                        " but it does not seem to exist!"))
            
    with open(filename, 'r') as f, open(temp_filename, 'w') as g :
        skip_next_lines = False
        for line in f :
            if line.startswith("F90 = ifort") :
                content = "F90 = ifort\n"\
                          "comp = -c -fPIC -O -module ${Mdir} -I${InDir}\n"\
                          "LFlagsB = -O -fPIC\n"\
                          "\n"\
                          "# Intels ifort, debug modus\n"\
                          "ifeq (${F90},ifortg)\n"\
                          " F90 = ifort\n"\
                          " comp = -c -g -fPIC -module ${Mdir} -I${InDir}\n"\
                          " LFlagsB = -g -fPIC\n"\
                          "endif\n"\
                          "\n"\
                          "# gfortran\n"\
                          "ifeq (${F90},gfortran)\n"\
                          " comp = -c -O -fPIC -J${Mdir} -I${InDir}\n"\
                          " LFlagsB = -O -fPIC\n"\
                          "endif\n"\
                          "\n"\
                          "# gfortran, any version (Added by GAMBIT)\n"\
                          "ifneq (,$(findstring gfortran,${F90}))\n"\
                          " comp = -c -O -fPIC -J${Mdir} -I${InDir}\n"\
                          " LFlagsB = -O -fPIC\n"\
                          "endif\n"\
                          "\n"\
                          "# g95 \n"\
                          "ifeq (${F90},g95)\n"\
                          " comp = -c -O -fPIC -fmod=${Mdir} -I${InDir}\n"\
                          " LFlagsB = -O -fPIC\n"\
                          "endif\n"\
                          "\n"\
                          "# Lahey F95 compiler\n"\
                          "ifeq (${F90},lf95)\n"\
                          " comp = -c -O -fPIC -M ${Mdir} -I${InDir}\n"\
                          " LFlagsB = -O -fPIC\n"\
                          "endif\n"\
                          "\n"\
                          "# NAG f95/2003\n"\
                          "ifeq (${F90},nagfor)\n"\
                          " comp = -c -O -fPIC -DONLYDOUBLE -mdir ${Mdir} -I${InDir}\n"\
                          " LFlagsB = -O -fPIC\n"\
                          "endif\n\n"
                g.write(content)
                skip_next_lines = True
            elif line.startswith(".SUFFIXES") :
                g.write(line)
                skip_next_lines = False
            else :
                if not skip_next_lines :
                    g.write(line)

            # We do not patch the -U flag in the ar command as it causes issues in OSX

    os.remove(filename)
    os.rename(temp_filename, filename)

def patch_control(model_name, patch_dir):
    """
    Patches $SPheno/src/Control.F90
    """
    
    filename = patch_dir + "/src/Control.F90"
    temp_filename = filename + "_temp"

    if not os.path.exists(filename):
        raise GumError(("Tried to find the file located at " + filename +
                        " but it does not seem to exist!"))
            
    with open(filename, 'r') as f, open(temp_filename, 'w') as g :
        skip_next_lines = False
        for line in f :
            if line.startswith(" Interface Is_NaN") :
                content = "! Added by GAMBIT\n"\
                          "Use, Intrinsic :: iso_c_binding\n"\
                          "Implicit none\n"\
                          "\n"\
                          "Type(c_funptr) :: ErrorHandler_cptr\n"\
                          "\n"\
                          "! Define interface of call-back routine.\n"\
                          "\n"\
                          "Abstract Interface\n"\
                          "  Subroutine callback ()\n"\
                          "    Use, Intrinsic :: iso_c_binding\n"\
                          "  End Subroutine callback\n"\
                          "End Interface\n"\
                          "\n"\
                          "! Variable to swith off screen output\n"\
                          "Logical :: SilenceOutput = .False.\n"\
                          "\n"\
                          "! GAMBIT addition end\n\n"
                g.write(content)
                g.write(line) 
            elif line.startswith(" Subroutine TerminateProgram") :
                content = " ! Subroutine modified by GAMBIT\n"\
                          " Subroutine TerminateProgram\n"\
                          " !-----------------------------------------------------------------------\n"\
                          " ! This subroutine terminates a program if a fatal error occurs.\n"\
                          " ! Before doing this, it writes the tree of calling subroutines to\n"\
                          " ! the file which is connected to the channel ErrCan\n"\
                          " ! written by Werner Porod, 20.9.2000\n"\
                          " !-----------------------------------------------------------------------\n"\
                          " Use, Intrinsic :: iso_c_binding\n"\
                          " Implicit None\n"\
                          " \n"\
                          "  Procedure(callback), Pointer :: ErrorHandler_fptr\n"\
                          "  Integer :: i1\n"\
                          " \n"\
                          "  Write (ErrCan,*) \"  \"\n"\
                          "  Write (ErrCan,*) \"ErrorLevel, Iname:\",ErrorLevel, Iname\n"\
                          "  Write (ErrCan,*) &\n"\
                          "    & \"The error has occured in the following chain of subroutines:\"\n"\
                          "  Do i1=1,Size(NameOfUnit)\n"\
                          "   Write (ErrCan,*) NameOfUnit(i1)\n"\
                          "  End Do\n"\
                          "  Write (ErrCan,*) \"  \"\n"\
                          "  Write (ErrCan,*) \"Hopefully you find the error soon\"\n"\
                          "  Do i1=ErrCan,ErrCan+NumberOfOpenFiles-1\n"\
                          "   Close(i1)\n"\
                          "  End Do\n"\
                          " \n"\
                          "  ! Convert C to Fortran procedure pointer.\n"\
                          "  Call c_f_procpointer(ErrorHandler_cptr, ErrorHandler_fptr)\n"\
                          " \n"\
                          "  ! Call the ErrorHandler\n"\
                          "  Call ErrorHandler_fptr()\n"\
                          " \n"\
                          "  ! This should never happen\n"\
                          "  Write(*,*) \"DEBUG: SPheno has continued past the ErrorHandler call. This should never happen...\"\n"\
                          "  Stop \"Subroutine TerminateProgram\"\n"\
                          " \n"\
                          " End Subroutine TerminateProgram\n\n"
                g.write(content)
                skip_next_lines = True
            elif line.startswith(" Subroutine") :
                g.write(line)
                skip_next_lines = False
            else :
                if not skip_next_lines :
                    g.write(line)


    os.remove(filename)
    os.rename(temp_filename, filename)


def patch_brs(model_name, patch_dir):
    """
    Patches $SPheno/<MODEL>/BranchingRatios_<MODEL>.f90
    """

    filename = patch_dir + "/" + model_name + "/BranchingRatios_" + model_name + ".f90"
    temp_filename = filename + "_temp"

    if not os.path.exists(filename):
        raise GumError(("Tried to find the file located at " + filename +
                        " but it does not seem to exist!"))
            
    with open(filename, 'r') as f, open(temp_filename, 'w') as g :
        for line in f :
            if line.startswith("Subroutine CalculateBR") :
                g.write(line[:22] + "_2" + line[22:])
            elif line.startswith("NameOfUnit(Iname) = \'CalculateBR\'") :
                content = "NameOfUnit(Iname) = \'CalculateBR_2\'\n"\
                                    "\n"\
                                    "! Added by GAMBIT\n"\
                                    "If (SilenceOutput) Then\n"\
                                    "  open(unit=6, file=\"/dev/null\", status=\"old\")\n"\
                                    "Endif\n"
                g.write(content)
            elif line.startswith("End Subroutine CalculateBR") :
                g.write("End Subroutine CalculateBR_2\n")
            else :
                g.write(line)


    os.remove(filename)
    os.rename(temp_filename, filename)


def patch_loopfunctions(model_name, patch_dir):
    """
    Patches $SPheno/<MODEL>/AddLoopFunctions.f90
    """

    filename = patch_dir + "/" + model_name + "/AddLoopFunctions.f90"
    temp_filename = filename + "_temp"

    if not os.path.exists(filename):
        raise GumError(("Tried to find the file located at " + filename +
                        " but it does not seem to exist!"))
            
    with open(filename, 'r') as f, open(temp_filename, 'w') as g :
        for line in f :
            if line.startswith("  SA_DerB00 = 3._dp * (xm1 * xm1 - 2._dp * xm1 * xm2") :
                content = "  If ((xm1.Eq.0._dp).And.(xm2.Eq.0._dp)) Then\n"\
                                    "    SA_DerB00 = -2._dp * xp * xp - 3._dp * xp * xp * VB0 - 3._dp * xp * xp * xp * VDerB0\n"\
                                    "  Else\n"
                g.write(content)
            elif line.startswith("  SA_DerB00 = SA_DerB00 / (36._dp * xp * xp)") :
                g.write("  Endif\n")
            g.write(line)

    os.remove(filename)
    os.rename(temp_filename, filename)


def patch_spheno_model(model_name, patch_dir):
    """
    Patches $SPheno/<MODEL>/SPheno<MODEL>.f90
    """
 
    filename = patch_dir + "/" + model_name + "/SPheno" + model_name + ".f90"
    temp_filename = filename + "_temp"

    if not os.path.exists(filename):
        raise GumError(("Tried to find the file located at " + filename +
                        " but it does not seem to exist!"))
            
    with open(filename, 'r') as f, open(temp_filename, 'w') as g:
        for line in f:
            if line.startswith("Program SPheno" + model_name) :
                g.write("!Program SPheno" + model_name + " ! Commented by GAMBIT\n")
                g.write("Module SPheno" + model_name + " ! Added by GAMBIT\n")
            elif line.startswith("Tpar = 0._dp") :
                content = "Contains ! Added by GAMBIT\n"\
                                    "\n"\
                                    "Subroutine Dummy() ! Added by GAMBIT\n"
                g.write(content)
                g.write(line)
            elif line.startswith(" Call CalculateBR") :
                g.write(line[:17] + "_2" + line[17:])
            elif line.startswith("Contains") :
                content = "\n"\
                                    "End Subroutine Dummy ! Added by GAMBIT\n"\
                                    "!Contains ! Commented by GAMBIT\n"
                g.write(content)
            elif line.startswith("kont = 0") :
                line2 = next(f)
                if line2.startswith("Call FirstGuess") :
                    content = "! Added by GAMBIT\n"\
                                        "If (SilenceOutput) Then\n"\
                                        " open(unit=6, file=\"/dev/null\", status=\"old\")\n"\
                                        "Endif\n"\
                                        "\n"\
                                        "kont = 0\n"
                    g.write(content)
                    g.write(line2)
            elif line.startswith("!If (kont.ne.0) Call TerminateProgram") :
                g.write("If (kont.ne.0) Call TerminateProgram\n")
            elif line.startswith("End Program SPheno" + model_name):
                g.write("!End Program SPheno" + model_name + " ! Commented by GAMBIT\n")
                g.write("End Module SPheno" + model_name + " ! Added by GAMBIT\n")
            else:
                g.write(line)

    os.remove(filename)
    os.rename(temp_filename, filename)


# SUSY-only patches

def patch_model_data(model_name, patch_dir):
    """
    Patches $SPheno/<MODEL>Model_Data_<MODEL>.f90
    """
    
    filename = patch_dir + "/" + model_name + "/Model_Data_" + model_name + ".f90"
    temp_filename = filename + "_temp"

    if not os.path.exists(filename):
        raise GumError(("Tried to find the file located at " + filename +
                        " but it does not seem to exist!"))
            
    with open(filename, 'r') as f, open(temp_filename, 'w') as g :
        for line in f :
            if line.startswith("Logical, Save :: CalcLoopDecay_LoopInducedOnly=.False.") :
                g.write(line)
                g.write("Logical, Save :: CalcSUSY3BodyDecays=.False. ! Added by GAMBIT\n")
            else :
                g.write(line)

    os.remove(filename)
    os.rename(temp_filename, filename)


def patch_3_body_decays_susy(model_name, patch_dir):
    """
    Patches the 3-body decays in: 
    $SPheno/<MODEL>/3-Body-Decays/X_<MODEL>.f90
    where X is a superfield.
    """

    particles = {"Cha", "Chi", "Glu", "Sd", "Su", "Se", "Sv"}
    channels = {"Cha": {"ChacChaCha", "ChaChiChi"},
                            "Chi": {"ChicChaCha", "ChiChiChi"},
                            "Glu": {},
                            "Sd": {"ChaGluSu", "SdChacCha", "SdChiChi", "ChiGluSd", "GluGluSd"},
                            "Su": {"SuChiChi", "ChiGluSu", "SdChicCha", "GluGluSu", "GluSdcCha", "SuChacCha"},
                            "Se": {"SvChaChi", "SeChacCha", "SeChiChi"},
                            "Sv": {"SvChiChi", "SeChicCha", "SvChacCha"}}
    # TODO: these channels work for the NMSSM, other susy models may have others

    for particle in particles :
 
        filename = patch_dir + "/" + model_name + "/3-Body-Decays/" + particle + "_" + model_name + ".f90"
        temp_filename = filename + "_temp"

        if not os.path.exists(filename):
            raise GumError(("Tried to find the file located at " + filename +
                            " but it does not seem to exist!"))
            
        with open(filename, 'r') as f, open(temp_filename, 'w') as g :
            for line in f :
                if line.startswith("Use ThreeBodyPhaseSpace") :
                    g.write(line)
                    g.write("Use Model_Data_" + model_name + " ! Added by GAMBIT\n")
                elif any([line.startswith("Call " + particle + "To" + channel) for channel in channels[particle]]) :
                    g.write("If (CalcSUSY3BodyDecays) Then ! Added by GAMBIT\n")
                    g.write(line)
                elif any([line.startswith("g" + particle + channel + "(i_run,:,:,:) = g" + particle + channel + "i") for channel in channels[particle]]) :
                    g.write("End If ! Added by GAMBIT\n\n")
                    g.write(line)
                else :
                    g.write(line)

        os.remove(filename)
        os.rename(temp_filename, filename)


"""
FRONTEND ROUTINES
"""


class SPhenoParameter:
    """
    Container type for a SPheno parameter.
    """
    
    def __init__(self, _name, _type, _size):

        self.name = _name
        self.type = _type
        self.size = _size


def write_spheno_frontends(model_name, parameters, spheno_path):
    """
    Writes the frontend source and header files for SPheno.
    """

    # Firstly, scrape the function signatures from the spheno source
    functions, arguments, locations = scrape_functions_from_spheno(spheno_path,
                                                                   model_name)

    # Convert the arguments to GAMBIT types
    type_dictionary = get_fortran_shapes(arguments)

    # Get all of the variables used in SPheno so we can store them as 
    # BE_VARIABLES. 
    variables = harvest_spheno_model_variables(spheno_path, model_name)

    # Get the source and header files
    spheno_src = write_spheno_frontend_src(model_name, functions)
    spheno_header = write_spheno_frontend_header(model_name, functions, 
                                                 type_dictionary, 
                                                 locations)


    #print spheno_header
    return spheno_src, spheno_header

def harvest_spheno_model_variables(spheno_path, model_name):
    """
    Harvests the model variables from $SPHENO/<MODEL>/Model_Data_<MODEL>.f90.
    Returns a dictionary of key: parameter name, value: GAMBIT fortran type, 
    and the same for the HiggsBounds output.
    """

    clean_model_name = model_name.replace('-','')
    location = "{0}/{1}/Model_Data_{1}.f90".format(spheno_path, 
                                                   clean_model_name)

        
    parameters = {}
    hb_parameters = {}

    src = ""
    hb_src = ""
    # First entry we care about is the line with "mass_uncertainty_Q" on it.
    # Last one will have the line "Contains". Read everything else in, then.
    # Save the HiggsBounds stuff too. This comes after a "For HiggsBounds" flag
    # and ends when the "HPPloop" variables are written, followed by 
    # " Real(dp) :: m32, tanbetaMZ "
    with open(location, 'r') as f:
        started = False
        hb_started = False
        for line in f:
            if "mass_uncertainty" in line:
                src += line
                started = True
            elif "HiggsBounds" in line:
                hb_started = True
            elif hb_started and "m32, tanbetaMZ" in line:
                hb_started = False
                started = False
            # Done.
            elif "Contains" in line:
                started = False
                hb_started = False
            elif started and not hb_started:
                src += line
            elif hb_started:
                hb_src += line

    # Source output -- clean it up a bit to make it easier to parse
    src = src.replace(' ','').replace('&\n',' ').replace('&','').split('\n')
    hb_src = hb_src.replace(' ','').replace('&\n',' ').replace('&','').split('\n')

    # The list of possible types a parameter could be
    possible_types = ["Real(dp)", "Integer", "Complex(dp)", "Logical"]

    # Each line looks like - TYPE :: definition(s)
    for line in src:
        # Each "split" is either an empty list, or a list with type first,
        # and parameter definition second.
        split = filter(None, line.split('::'))
        if not split: continue # Empty list -- skip it
        # Get the type.
        _type = ""
        for pt in possible_types:
            if pt in split[0]:
                _type = pt
        if _type == "":
            raise GumError("No type scraped from Model_Data.")

        # Split the RHS up but *not* if there is a comma between parentheses.
        names = re.split(r',\s*(?![^()]*\))', split[1])
        for name in names:
            # If it's defaulted, don't wanna save that.
            if "=" in name:
                name = name.split('=')[0]
            # Any size information? Like (3,3) or (2)...
            pat = '\((.*?)\)'
            r = re.search(pat, name)
            if r:
                size = r.group(1)
                name = name.split('(')[0] # Save the name without the (..)
            else: 
                size = ""

            par = SPhenoParameter(name, _type, size)
            parameters[name] = par

    # Exactly the same for HiggsBounds output
    for line in hb_src:
        split = filter(None, line.split('::'))
        if not split: continue 
        _type = ""
        for pt in possible_types:
            if pt in split[0]:
                _type = pt
        if _type == "":
            raise GumError("No type scraped from Model_Data.")
        names = re.split(r',\s*(?![^()]*\))', split[1])
        for name in names:
            if "=" in name:
                name = name.split('=')[0]
            pat = '\((.*?)\)'
            r = re.search(pat, name)
            if r:
                size = r.group(1)
                name = name.split('(')[0] # Save the name without the (..)
            else: 
                size = ""

            par = SPhenoParameter(name, _type, size)
            hb_parameters[name] = par
    
    return parameters, hb_parameters

def get_fortran_shapes(parameters):
    """
    Returns the GAMBIT fortran-shape for a SPheno parameter
    """

    type_dictionary = {}

    for name, parameter in parameters.iteritems():

        if not isinstance(parameter, SPhenoParameter):
            raise GumError(("Parameter not passed as instance of "
                            "SPhenoParameter to function get_fortran_shape."))

        fortran_type = ""
        
        # If it is not a scalar, then it's an Farray
        if parameter.size:
            fortran_type += "Farray_"

        # If it's complex
        if parameter.type == "Complex(dp)":
            fortran_type += "Fcomplex16"

        elif parameter.type == "Real(dp)":
            fortran_type += "Freal8"

        elif parameter.type == "Integer":
            fortran_type += "Finteger"

        elif parameter.type == "Logical":
            fortran_type += "Flogical"

        # Now define the size of the array, if there is one
        if parameter.size:
            fortran_type += "_1_"
            fortran_type += parameter.size.replace(',','_1_')

        # If we haven't been able to convert it, throw an error.
        if not fortran_type:
            raise GumError(("GUM has not been able to convert the type of the "
                            "parameter " + name + ". Please check "
                            "your SPheno file, and if necessary, add the type "
                            "to get_fortran_shape."))

        type_dictionary[name] = fortran_type  

    return type_dictionary



def get_arguments_from_file(functions, file_path, function_dictionary,
                            argument_dictionary):
    """
    Helper function to obtain the function signatues from a given file.
    Also pull the types of each argument while we are there.
    """

    # The list of possible types a parameter could be
    possible_types = ["Real(dp)", "Integer", "Complex(dp)", "Logical"]

    with open(file_path) as f:
        data = f.readlines()
    # Remove line breaks -- makes regex easier to use
    file = "".join(data).replace("\n","")

    # If you pass a string instead of a list, wrap it into one
    if not isinstance(functions, list):
        funcs = []
        funcs.append(functions)
        functions = funcs

    for function in functions:

        # Extract the signature for a given function
        pattern = 'Subroutine {0}\((.*?)\)'.format(function)
        s = re.search(pattern, file)
        if s:
            signature = s.group(1).replace('&','').replace(' ','').split(',')
            function_dictionary[function] = signature
        else:
            raise GumError("Error extracting " + function +" from SPheno.")

        # Get some function definitions while we are here too.

        # This is the entire function definition we're grabbing
        pattern = 'Subroutine {0}(.*?)End Subroutine {0}'.format(function)
        s = re.search(pattern, file)
        if s:
            func = s.group(1)
            # Split the definition up by the :: sign for definitions.
            defs = func.split('::')



            # Firstly, find out which arguments we need to obtain types for
            for v in function_dictionary.values():
                for arg in v:
                    # If we've got it already
                    if arg in argument_dictionary: 
                        continue
                    else:
                        # Go through each split; if we find a possible Fortran type
                        # in the string, scrape the parameter name and add to the 
                        # dictionary.
                        for i in range(len(defs)):
                            # ... unless we've already got it.
                            if arg in argument_dictionary: 
                                continue
                            else:
                                for j in possible_types:
                                    if j in defs[i]:
                                        if arg in defs[i+1]:
                                            # Also check the size here by looking
                                            # at the size of the matrix
                                            pat = '{}\((.*?)\)'.format(arg)
                                            r = re.search(pat, defs[i+1])
                                            if r:
                                                size = r.group(1)
                                            else: 
                                                size = ""
                                            # Size can also be given by "Dimension"
                                            # in FORTRAN
                                            if "Dimension" in defs[i]:
                                                size = re.search(r'Dimension\((.*?)\)', defs[i]).group(1)
                                            par = SPhenoParameter(arg, j, size)
                                            argument_dictionary[arg] = par
                                            break

        else:
            raise GumError(("Couldn't match the pattern between the definition"
                            " and undefinition of the Subroutine " + function +
                            ", please inspect the SARAH-SPheno source code, "
                            "as I don't think it will work anyway."))



def scrape_functions_from_spheno(spheno_path, model_name):
    """
    Reads the SPheno source to identify the function signatures
    we want to include in the frontend.
    """
   
    clean_model_name = model_name.replace('-','')

    # Create a dictionary of function names, and the string of 
    # parameters that go into them.
    function_dictionary = {}
    args_dictionary = {}

    # Dictionary of which files each function lives in.
    locations_dictionary = {
        clean_model_name+"/BranchingRatios_"+clean_model_name  : "CalculateBR_2",
        clean_model_name+"/Unitarity_"+clean_model_name        : ["ScatteringEigenvalues",
                                                                  "ScatteringEigenvaluesWithTrilinears"],
        clean_model_name+"/SPheno"+clean_model_name            : ["CalculateSpectrum",
                                                                  "GetScaleUncertainty"],
        clean_model_name+"/Model_Data_"+clean_model_name       : "SetMatchingConditions",
        clean_model_name+"/LoopMasses_"+clean_model_name       : "OneLoopMasses",
        clean_model_name+"/InputOutput_"+clean_model_name      : ["Switch_to_superCKM",
                                                                  "Switch_to_superPMNS"],
        "src/Model_Data"                                       : ["Switch_from_superCKM",
                                                                  "Switch_from_superPMNS"],
        clean_model_name+"/TadpoleEquations_"+clean_model_name : "SolveTadpoleEquations"
    }

    for location, functions in locations_dictionary.iteritems():
        filename = spheno_path+"/"+location+".f90"
        get_arguments_from_file(functions, filename, function_dictionary, args_dictionary)

    return function_dictionary, args_dictionary, locations_dictionary

def write_spheno_frontend_src(model_name, function_signatures):
    """
    Writes source for 
    Backends/src/frontends/SARAHSPheno_<MODEL>_<VERSION>.cpp
    """ 
    intro_message = "Frontend for SARAH-SPheno {0} backend, for the "\
                    " {1} model.".format(SPHENO_VERSION, model_name)

    towrite = blame_gum(intro_message)

    towrite += (
        "#include \"gambit/Backends/frontend_macros.hpp\"\n"
        "#include \"gambit/Backends/frontends/SARAHSPheno_{0}_{1}.hpp\"\n"
        "#include \"gambit/Elements/spectrum_factories.hpp\"\n"
        "#include \"gambit/Models/SimpleSpectra/{0}SimpleSpec.hpp\"\n"
        "#include \"gambit/Utils/version.hpp\"\n"
        "\n"
        "#define BACKEND_DEBUG 0\n"
        "\n"
        "// Callback function for error handling\n"
        "BE_NAMESPACE\n"
        "{{\n"
        "// This function will be called from SPheno. Needs C linkage, and thus also\n"
        "// a backend-specific name to guard against name clashes.\n"
        "extern \"C\"\n"
        "void CAT_4(BACKENDNAME,_,SAFE_VERSION,_ErrorHandler)()\n"
        "{{\n"
        "throw std::runtime_error(\"SARAHSPheno_{0} backend called TerminateProgram.\");\n"
        "}}\n"
        "}}\n"
        "END_BE_NAMESPACE\n"
        "\n"
        "// Convenience functions (definition)\n"
        "BE_NAMESPACE\n"
        "{{\n"
        "\n"
        "// Variables and functions to keep and access decay info\n"
        "typedef std::tuple<std::vector<int>,int,double> channel_info_triplet; // {{pdgs of daughter particles}}, spheno index, correction factor\n"
        "namespace Fdecays\n"
        "{{\n"
        "// A (pdg,vector) map, where the vector contains a channel_info_triplet for each\n"
        "// decay mode of the mother particle. (See typedef of channel_info_triplet above.)\n"
        "static std::map<int,std::vector<channel_info_triplet> > all_channel_info;\n"
        "\n"
        "// Flag indicating whether the decays need to be computed or not.\n"
        "static bool BRs_already_calculated = false;\n"
        "\n"
        "// Function that reads a table of all the possible decays in SARAHSPheno_{0}\n"
        "// and fills the all_channel_info map above\n"
        "void fill_all_channel_info(str);\n"
        "\n"
        "// Helper function to turn a vector<int> into a vector<pairs<int,int> > needed for\n"
        "// when calling the GAMBIT DecayTable::set_BF function\n"
        "std::vector<std::pair<int,int> > get_pdg_context_pairs(std::vector<int>);\n"
        "}}\n"
        "\n"
    ).format(model_name, SPHENO_VERSION.replace('.','_'))


    #print indent(towrite)
    return indent(towrite)

def write_spheno_frontend_header(model_name, function_signatures, type_dictionary, locations):
    """
    Writes code for 
    Backends/include/gambit/Backends/SARAHSPheno_<MODEL>_<VERSION>.hpp

    Information needed:
        -> parameters (as they are known in SPheno) and their type and size,
        and the order they appear in the functions.
        -> function signatures from spheno
        -> number of Higgses
        -> if it's a SUSY model
    """

    clean_model_name = model_name.replace('-','')

    intro_message = "Frontend header for SARAH-SPheno {0} backend, for the "\
                    " {1} model.".format(SPHENO_VERSION, model_name)
                    
    towrite = blame_gum(intro_message)

    # Some nice model-independent functions to begin
    towrite += (
            "#define BACKENDNAME SARAHSPheno_{0}\n"
            "#define BACKENDLANG FORTRAN\n"
            "#define VERSION {1}\n"
            "#define SARAH_VERSION {2}\n"
            "#define SAFE_VERSION {3}\n"
            "\n"
            "// Begin\n"
            "LOAD_LIBRARY\n"
            "\n"
            "// Allow for {0} only\n"
            "BE_ALLOW_MODELS({0})\n"
            "\n"
            "\n"
            "// Functions\n"
            "BE_FUNCTION(Set_All_Parameters_0, void, (), \"__model_data_{4}_MOD_set_all_parameters_0\", \"SARAHSPheno_{0}_internal\")\n"
            "BE_FUNCTION(SetRenormalizationScale, Freal8, (Freal8&), \"__loopfunctions_MOD_setrenormalizationscale\", \"SARAHSPheno_{0}_internal\")\n"
            "BE_FUNCTION(InitializeLoopFunctions, void, (), \"__loopfunctions_MOD_initializeloopfunctions\", \"SARAHSPheno_{0}_internal\")\n"
            "BE_FUNCTION(CalculateRunningMasses, void, (Farray_Freal8_1_3&, //mf_l_in\n"
            "                                           Farray_Freal8_1_3&, // mf_d_in\n"
            "                                           Farray_Freal8_1_3&, // mf_u_in\n"
            "                                           Freal8&, // Qlow\n"
            "                                           Freal8&, // Alpha\n"
            "                                           Freal8&, // AlphaS\n"
            "                                           Freal8&, // Qhigh\n"
            "                                           Farray_Freal8_1_3&, // mf_l_out\n"
            "                                           Farray_Freal8_1_3&, // mf_d_out\n"
            "                                           Farray_Freal8_1_3&, // mf_u_out\n"
            "                                           Finteger&), //kont))\n"
            "     \"__standardmodel_MOD_calculaterunningmasses\", \"SARAHSPheno_{0}_internal\")\n"
            "BE_FUNCTION(GetRenormalizationScale, Freal8, (), \"__loopfunctions_MOD_getrenormalizationscale\", \"SARAHSPheno_{0}_internal\")\n"
            "BE_FUNCTION(SetRGEScale, void, (Freal8&), \"__model_data_{4}_MOD_setrgescale\", \"SARAHSPheno_{0}_internal\")\n"
            "BE_FUNCTION(SetGUTScale, void, (Freal8&), \"__model_data_{4}_MOD_setgutscale\", \"SARAHSPheno_{0}_internal\")\n"
            "BE_FUNCTION(SetStrictUnification, Flogical, (Flogical&), \"__model_data_{4}_MOD_setstrictunification\", \"SARAHSPheno_{0}_internal\")\n"
            "BE_FUNCTION(SetYukawaScheme, Finteger, (Finteger&), \"__model_data_{4}_MOD_setyukawascheme\", \"SARAHSPheno_{0}_internal\")\n"
    ).format(clean_model_name, SPHENO_VERSION, SARAH_VERSION, SPHENO_VERSION.replace('.','_'),
             clean_model_name.lower())

    # Some model-dependent functions:
    """
    ScatteringEigenvalues; GetScaleUncertainty; CalculateSpectrum; OneLoopMasses;
    SolveTadpoleEquations; Switch_to_superPMNS; Switch_to_superCKM; Switch_from_superCKM;
    SetMatchingConditions; ScatteringEigenvaluesWithTrilinears; CalculateBR
    """

    towrite += "\n// Model-dependent arguments auto-scraped by GUM\n"
    
    # The functions that we have scraped from SPheno directly.
    for function, sig in function_signatures.iteritems():
        # And the locations where each function lives
        for k, v in locations.iteritems():
            # If they match, store the module name (usually the filename)
            # as this is going to be the name of the symbol in SPheno
            if function in v:
                loc = k.lower().split('/')[-1]
                # Overwrite the tadpoles module, the filename is different
                if "tadpole" in loc: loc = "tadpoles"
                symbol = "__{0}_MOD_{2}".format(loc, 
                         clean_model_name.lower(), function.lower())
        # Now list all arguments, with a nice comment next to it, to make it 
        # lovely and readable.
        arguments = []
        arguments.append(type_dictionary[sig[0]] + '&, // ' + sig[0])
        for argument in sig[1:-1]:
            arguments.append("   "+type_dictionary[argument]+"&, // "+argument)
        arguments.append("   "+type_dictionary[sig[-1]]+"& // "+sig[-1]+"\n")
        args = "\n".join(arguments)
        # Wrap it all up.
        towrite += (
            "BE_FUNCTION({0}, void,\n"
            "  ({1}),"
            " \"{2}\", \"SARAHSPheno_{3}_internal\")\n"
        ).format(function, args, symbol, clean_model_name)
    
    # MODEL VARIABLES
    # TODO

    # SPHENOINPUT VARIABLES
    # TODO

    # MINPAR VARIABLES
    # TODO
    
    # EXTPAR VARIABLES
    # TODO

    # SMINPUTS
    towrite += (
            "// SMINPUT Variables\n"
            "BE_VARIABLE(mZ, Freal8, \"__standardmodel_MOD_mz\", \"SARAHSPheno_{0}_internal\")\n"
            "BE_VARIABLE(mZ2, Freal8,  \"__standardmodel_MOD_mz2\", \"SARAHSPheno_{0}_internal\")\n"
            "BE_VARIABLE(gamZ, Freal8, \"__standardmodel_MOD_gamz\", \"SARAHSPheno_{0}_internal\")\n"
            "BE_VARIABLE(gamZ2, Freal8, \"__standardmodel_MOD_gamz2\", \"SARAHSPheno_{0}_internal\")\n"
            "BE_VARIABLE(gmZ, Freal8, \"__standardmodel_MOD_gmz\", \"SARAHSPheno_{0}_internal\")\n"
            "BE_VARIABLE(gmZ2, Freal8, \"__standardmodel_MOD_gmz2\", \"SARAHSPheno_{0}_internal\")\n"
            "BE_VARIABLE(BrZqq, Farray_Freal8_1_5, \"__standardmodel_MOD_brzqq\", \"SARAHSPheno_{0}_internal\")\n"
            "BE_VARIABLE(BrZll, Farray_Freal8_1_3, \"__standardmodel_MOD_brzll\", \"SARAHSPheno_{0}_internal\")\n"
            "BE_VARIABLE(BrZinv, Freal8, \"__standardmodel_MOD_brzinv\", \"SARAHSPheno_{0}_internal\")\n"
            "BE_VARIABLE(mW, Freal8, \"__standardmodel_MOD_mw\", \"SARAHSPheno_{0}_internal\")\n"
            "BE_VARIABLE(mW_SM, Freal8, \"__model_data_{0}_MOD_mw_sm\", \"SARAHSPheno_{0}_internal\")\n"
            "BE_VARIABLE(mW2, Freal8, \"__standardmodel_MOD_mw2\", \"SARAHSPheno_{0}_internal\")\n"
            "BE_VARIABLE(gamW, Freal8, \"__standardmodel_MOD_gamw\", \"SARAHSPheno_{0}_internal\")\n"
            "BE_VARIABLE(gamW2, Freal8, \"__standardmodel_MOD_gamw2\", \"SARAHSPheno_{0}_internal\")\n"
            "BE_VARIABLE(gmW, Freal8, \"__standardmodel_MOD_gmw\", \"SARAHSPheno_{0}_internal\")\n"
            "BE_VARIABLE(gmW2, Freal8, \"__standardmodel_MOD_gmw2\", \"SARAHSPheno_{0}_internal\")\n"
            "BE_VARIABLE(BrWqq, Farray_Freal8_1_2, \"__standardmodel_MOD_brwqq\", \"SARAHSPheno_{0}_internal\")\n"
            "BE_VARIABLE(BrWln, Farray_Freal8_1_3, \"__standardmodel_MOD_brwln\", \"SARAHSPheno_{0}_internal\")\n"
            "BE_VARIABLE(mf_l, Farray_Freal8_1_3, \"__standardmodel_MOD_mf_l\", \"SARAHSPheno_{0}_internal\")\n"
            "BE_VARIABLE(mf_l_mZ, Farray_Freal8_1_3, \"__standardmodel_MOD_mf_l_mz\", \"SARAHSPheno_{0}_internal\")\n"
            "BE_VARIABLE(mf_nu, Farray_Freal8_1_3, \"__standardmodel_MOD_mf_nu\", \"SARAHSPheno_{0}_internal\")\n"
            "BE_VARIABLE(mf_u, Farray_Freal8_1_3, \"__standardmodel_MOD_mf_u\", \"SARAHSPheno_{0}_internal\")\n"
            "BE_VARIABLE(mf_u_mZ, Farray_Freal8_1_3, \"__standardmodel_MOD_mf_u_mz\", \"SARAHSPheno_{0}_internal\")\n"
            "BE_VARIABLE(mf_d, Farray_Freal8_1_3, \"__standardmodel_MOD_mf_d\", \"SARAHSPheno_{0}_internal\")\n"
            "BE_VARIABLE(mf_d_mZ, Farray_Freal8_1_3, \"__standardmodel_MOD_mf_d_mz\", \"SARAHSPheno_{0}_internal\")\n"
            "BE_VARIABLE(mf_l2, Farray_Freal8_1_3, \"__standardmodel_MOD_mf_l2\", \"SARAHSPheno_{0}_internal\")\n"
            "BE_VARIABLE(mf_u2, Farray_Freal8_1_3, \"__standardmodel_MOD_mf_u2\", \"SARAHSPheno_{0}_internal\")\n"
            "BE_VARIABLE(mf_d2, Farray_Freal8_1_3, \"__standardmodel_MOD_mf_d2\", \"SARAHSPheno_{0}_internal\")\n"
            "BE_VARIABLE(MNuR, Freal8, \"__model_data_MOD_mnur\", \"SARAHSPheno_{0}_internal\")\n"
            "BE_VARIABLE(Q_light_quarks, Freal8, \"__standardmodel_MOD_q_light_quarks\", \"SARAHSPheno_{0}_internal\")\n"
            "BE_VARIABLE(Delta_Alpha_Lepton, Freal8, \"__standardmodel_MOD_delta_alpha_lepton\", \"SARAHSPheno_{0}_internal\")\n"
            "BE_VARIABLE(Delta_Alpha_Hadron, Freal8, \"__standardmodel_MOD_delta_alpha_hadron\", \"SARAHSPheno_{0}_internal\")\n"
            "BE_VARIABLE(Alpha, Freal8, \"__standardmodel_MOD_alpha\", \"SARAHSPheno_{0}_internal\")\n"
            "BE_VARIABLE(Alpha_mZ, Freal8, \"__standardmodel_MOD_alpha_mz\", \"SARAHSPheno_{0}_internal\")\n"
            "BE_VARIABLE(Alpha_mZ_MS, Freal8, \"__standardmodel_MOD_alpha_mz_ms\", \"SARAHSPheno_{0}_internal\")\n"
            "BE_VARIABLE(MZ_input, Flogical, \"__model_data_{0}_MOD_mz_input\", \"SARAHSPheno_{0}_internal\")\n"
            "BE_VARIABLE(AlphaS_mZ, Freal8, \"__standardmodel_MOD_alphas_mz\", \"SARAHSPheno_{0}_internal\")\n"
            "BE_VARIABLE(G_F, Freal8, \"__standardmodel_MOD_g_f\", \"SARAHSPheno_{0}_internal\")\n"
            "BE_VARIABLE(KFactorLee, Freal8, \"__standardmodel_MOD_kfactorlee\", \"SARAHSPheno_{0}_internal\")\n"
            "BE_VARIABLE(CKM, Farray_Fcomplex16_1_3_1_3, \"__standardmodel_MOD_ckm\", \"SARAHSPheno_{0}_internal\")\n"
            "BE_VARIABLE(lam_wolf, Freal8, \"__standardmodel_MOD_lam_wolf\", \"SARAHSPheno_{0}_internal\")\n"
            "BE_VARIABLE(A_wolf, Freal8, \"__standardmodel_MOD_a_wolf\", \"SARAHSPheno_{0}_internal\")\n"
            "BE_VARIABLE(rho_wolf, Freal8, \"__standardmodel_MOD_rho_wolf\", \"SARAHSPheno_{0}_internal\")\n"
            "BE_VARIABLE(eta_wolf, Freal8, \"__standardmodel_MOD_eta_wolf\", \"SARAHSPheno_{0}_internal\")\n"
            "\n"
    ).format(clean_model_name)

    # MASS + OUTPUT VARIABLES
    # TODO

    # MODEL VARIABLES
    # TODO

    # CONTROL  + "OTHER" VARIABLES
    towrite += (
            "// Control Variables\n"
            "BE_VARIABLE(Iname, Finteger, \"__control_MOD_iname\", \"SARAHSPheno_{0}_internal\")\n"
            "BE_VARIABLE(kont, Finteger, \"__spheno{1}_MOD_kont\", \"SARAHSPheno_{0}_internal\")\n"
            "BE_VARIABLE(WriteOut, Flogical, \"__control_MOD_writeout\", \"SARAHSPheno_{0}_internal\")\n"
            "BE_VARIABLE(epsI, Freal8, \"__spheno{1}_MOD_epsi\", \"SARAHSPheno_{0}_internal\")\n"
            "BE_VARIABLE(deltaM, Freal8, \"__spheno{1}_MOD_deltam\", \"SARAHSPheno_{0}_internal\")\n"
            "BE_VARIABLE(mGUT, Freal8, \"__spheno{1}_MOD_mgut\", \"SARAHSPheno_{0}_internal\")\n"
            "BE_VARIABLE(ErrCan, Finteger, \"__control_MOD_errcan\", \"SARAHSPheno_{0}_internal\")\n"
            "BE_VARIABLE(FoundIterativeSolution, Flogical, \"__settings_MOD_founditerativesolution\", \"SARAHSPheno_{0}_internal\")\n"
            "\n"   
            "// Other variables\n"
            "BE_VARIABLE(Qin, Freal8, \"__spheno{1}_MOD_qin\", \"SARAHSPheno_{0}_internal\")\n"
            "BE_VARIABLE(ratioWoM, Freal8, \"__spheno{1}_MOD_ratiowom\",\"SARAHSPheno_{0}_internal\")\n"
            "BE_VARIABLE(CalcTBD, Flogical, \"__spheno{1}_MOD_calctbd\",\"SARAHSPheno_{0}_internal\")\n"
    ).format(clean_model_name, clean_model_name.lower())

    # BRANCHING RATIOS
    towrite += (
            "// Branching Ratio variables\n"
            "BE_VARIABLE(L_BR, Flogical, \"__control_MOD_l_br\", \"SARAHSPheno_{0}_internal\")\n"
            "BE_VARIABLE(Enable3BDecaysF, Flogical, \"__settings_MOD_enable3bdecaysf\", \"SARAHSPheno_{0}_internal\")\n"
            "BE_VARIABLE(Enable3BDecaysS, Flogical, \"__settings_MOD_enable3bdecayss\", \"SARAHSPheno_{0}_internal\")\n"
            "BE_VARIABLE(RunningCouplingsDecays, Flogical, \"__settings_MOD_runningcouplingsdecays\", \"SARAHSPheno_{0}_internal\")\n"
            "BE_VARIABLE(MinWidth, Freal8, \"__settings_MOD_minwidth\", \"SARAHSPheno_{0}_internal\")\n"
            "BE_VARIABLE(OneLoopDecays, Flogical, \"__settings_MOD_oneloopdecays\", \"SARAHSPheno_{0}_internal\")\n"
    ).format(clean_model_name)

    # DECAY OPTIONS
    towrite += (
            "\n"
            "// Decay options"
            "BE_VARIABLE(divonly_save, Finteger, \"__settings_MOD_divonly_save\", \"SARAHSPheno_{0}_internal\")\n"
            "BE_VARIABLE(divergence_save, Freal8, \"__settings_MOD_divergence_save\", \"SARAHSPheno_{0}_internal\")\n"
            "BE_VARIABLE(SimplisticLoopDecays, Flogical, \"__settings_MOD_simplisticloopdecays\", \"SARAHSPheno_{0}_internal\")\n"
            "BE_VARIABLE(ShiftIRdiv, Flogical, \"__settings_MOD_shiftirdiv\", \"SARAHSPheno_{0}_internal\")\n"
            "BE_VARIABLE(DebugLoopDecays, Flogical, \"__settings_MOD_debugloopdecays\", \"SARAHSPheno_{0}_internal\")\n"
            "BE_VARIABLE(OnlyTreeLevelContributions, Flogical, \"__settings_MOD_onlytreelevelcontributions\", \"SARAHSPheno_{0}_internal\")\n"
            "BE_VARIABLE(ExternalZfactors, Flogical, \"__settings_MOD_externalzfactors\", \"SARAHSPheno_{0}_internal\")\n"
            "BE_VARIABLE(UseZeroRotationMatrices, Flogical, \"__settings_MOD_usezerorotationmatrices\", \"SARAHSPheno_{0}_internal\")\n"
            "BE_VARIABLE(UseP2Matrices, Flogical, \"__settings_MOD_usep2matrices\", \"SARAHSPheno_{0}_internal\")\n"
            "BE_VARIABLE(OSkinematics, Flogical, \"__settings_MOD_oskinematics\", \"SARAHSPheno_{0}_internal\")\n"
            "BE_VARIABLE(ewOSinDecays, Flogical, \"__settings_MOD_ewosindecays\", \"SARAHSPheno_{0}_internal\")\n"
            "BE_VARIABLE(yukOSinDecays, Flogical, \"__settings_MOD_yukosindecays\", \"SARAHSPheno_{0}_internal\")\n"
            "BE_VARIABLE(CTinLoopDecays, Flogical, \"__settings_MOD_ctinloopdecays\", \"SARAHSPheno_{0}_internal\")\n"
            "BE_VARIABLE(LoopInducedDecaysOS, Flogical, \"__settings_MOD_loopinduceddecaysos\", \"SARAHSPheno_{0}_internal\")\n"
            "BE_VARIABLE(Mass_Regulator_PhotonGluon, Freal8, \"__settings_MOD_mass_regulator_photongluon\", \"SARAHSPheno_{0}_internal\")\n"
            "BE_VARIABLE(Extra_Scale_LoopDecays, Flogical, \"__settings_MOD_extra_scale_loopdecays\", \"SARAHSPheno_{0}_internal\")\n"
            "BE_VARIABLE(Scale_LoopDecays, Freal8, \"__settings_MOD_scale_loopdecays\", \"SARAHSPheno_{0}_internal\")\n"
    ).format(clean_model_name)

    # TODO: if is_susy:
    if "MSSM" in model_name: 
        towrite += (
                "BE_VARIABLE(Calc3BodyDecay_Glu, Flogical, \"__model_data_{0}_MOD_calc3bodydecay_glu\", \"SARAHSPheno_{1}_internal\")\n"
                "BE_VARIABLE(Calc3BodyDecay_Chi, Flogical, \"__model_data_{0}_MOD_calc3bodydecay_chi\", \"SARAHSPheno_{1}_internal\")\n"
                "BE_VARIABLE(Calc3BodyDecay_Cha, Flogical, \"__model_data_{0}_MOD_calc3bodydecay_cha\", \"SARAHSPheno_{1}_internal\")\n"
                "BE_VARIABLE(Calc3BodyDecay_Sd, Flogical, \"__model_data_{0}_MOD_calc3bodydecay_sd\", \"SARAHSPheno_{1}_internal\")\n"
                "BE_VARIABLE(Calc3BodyDecay_Su, Flogical, \"__model_data_{0}_MOD_calc3bodydecay_su\", \"SARAHSPheno_{1}_internal\")\n"
                "BE_VARIABLE(Calc3BodyDecay_Se, Flogical, \"__model_data_{0}_MOD_calc3bodydecay_se\", \"SARAHSPheno_{1}_internal\")\n"
                "BE_VARIABLE(Calc3BodyDecay_Sv, Flogical, \"__model_data_{0}_MOD_calc3bodydecay_sv\", \"SARAHSPheno_{1}_internal\")\n"
                "BE_VARIABLE(CalcSUSY3BodyDecays, Flogical, \"__model_data_{0}_MOD_calcsusy3bodydecays\", \"SARAHSPheno_{1}_internal\")\n"
                "BE_VARIABLE(CalcLoopDecay_LoopInducedOnly, Flogical, \"__model_data_{0}_MOD_calcloopdecay_loopinducedonly\", \"SARAHSPheno_{1}_internal\")\n"
                "BE_VARIABLE(CalcLoopDecay_Sd, Flogical, \"__model_data_{0}_MOD_calcloopdecay_sd\", \"SARAHSPheno_{1}_internal\")\n"
                "BE_VARIABLE(CalcLoopDecay_Su, Flogical, \"__model_data_{0}_MOD_calcloopdecay_su\", \"SARAHSPheno_{1}_internal\")\n"
                "BE_VARIABLE(CalcLoopDecay_Se, Flogical, \"__model_data_{0}_MOD_calcloopdecay_se\", \"SARAHSPheno_{1}_internal\")\n"
                "BE_VARIABLE(CalcLoopDecay_Sv, Flogical, \"__model_data_{0}_MOD_calcloopdecay_sv\", \"SARAHSPheno_{1}_internal\")\n"
                "BE_VARIABLE(CalcLoopDecay_hh, Flogical, \"__model_data_{0}_MOD_calcloopdecay_hh\", \"SARAHSPheno_{1}_internal\")\n"
                "BE_VARIABLE(CalcLoopDecay_Ah, Flogical, \"__model_data_{0}_MOD_calcloopdecay_ah\", \"SARAHSPheno_{1}_internal\")\n"
                "BE_VARIABLE(CalcLoopDecay_Hpm, Flogical, \"__model_data_{0}_MOD_calcloopdecay_hpm\", \"SARAHSPheno_{1}_internal\")\n"
                "BE_VARIABLE(CalcLoopDecay_Glu, Flogical, \"__model_data_{0}_MOD_calcloopdecay_glu\", \"SARAHSPheno_{1}_internal\")\n"
                "BE_VARIABLE(CalcLoopDecay_Chi, Flogical, \"__model_data_{0}_MOD_calcloopdecay_chi\", \"SARAHSPheno_{1}_internal\")\n"
                "BE_VARIABLE(CalcLoopDecay_Cha, Flogical, \"__model_data_{0}_MOD_calcloopdecay_cha\", \"SARAHSPheno_{1}_internal\")\n"
                "BE_VARIABLE(CalcLoopDecay_Fu, Flogical, \"__model_data_{0}_MOD_calcloopdecay_fu\", \"SARAHSPheno_{1}_internal\")\n"
        ).format(clean_model_name.lower(), clean_model_name)


    # HIGGSBOUNDS OUTPUT
    # TODO

    # Wrap it up.
    towrite += (
            "// Convenience functions (registration)\n"
            "BE_CONV_FUNCTION(run_SPheno, int, (Spectrum&, const Finputs&), \"{0}_spectrum\")\n"
            "BE_CONV_FUNCTION(run_SPheno_decays, int, (const Spectrum &, DecayTable &, const Finputs&), \"{0}_decays\")\n"
            "BE_CONV_FUNCTION(Spectrum_Out, Spectrum, (const std::map<str, safe_ptr<double> >&), \"SARAHSPheno_{0}_internal\")\n"
            "BE_CONV_FUNCTION(get_HiggsCouplingsTable, int, (const Spectrum&, HiggsCouplingsTable&, const Finputs&), \"{0}_HiggsCouplingsTable\")\n"
            "BE_CONV_FUNCTION(ReadingData, void, (const Finputs&), \"SARAHSPheno_{0}_internal\")\n"
            "BE_CONV_FUNCTION(ReadingData_decays, void, (const Finputs&), \"SARAHSPheno_{0}_internal\")\n"
            "BE_CONV_FUNCTION(InitializeStandardModel, void, (const SMInputs&), \"SARAHSPheno_{0}_internal\")\n"
            "BE_CONV_FUNCTION(ErrorHandling, void, (const int&), \"SARAHSPheno_{0}_internal\")\n"
            "\n"
            "// Initialisation functions (dependencies)\n"
            "\n"
            "// Function pointer variable for error handling\n"
            "BE_VARIABLE(ErrorHandler_cptr, type_fptr_SPhenoErrorHandler, \"__control_MOD_errorhandler_cptr\", \"SARAHSPheno_{0}_internal\")\n"
            "\n"
            "// End\n"
            "#include \"gambit/Backends/backend_undefs.hpp\"\n"
    ).format(clean_model_name)

    #print indent(towrite)
    return indent(towrite)


"""
SPECBIT ROUTINES
"""

"""
DECAYBIT ROUTINES
"""