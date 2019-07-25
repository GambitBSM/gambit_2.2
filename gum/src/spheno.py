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
from distutils.dir_util import copy_tree
from collections import defaultdict
import re

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

  # TODO: throw some errors if the file does not exist

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

  # TODO: throw some errors if the file does not exist

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
                  "\tmv SPhenoNMSSM ../bin\n"\
                  "\trm SPhenoNMSSM.o\n"
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

  # TODO: throw some errors if the file does not exist

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

      # We do not patch the -U flag in the ar command as it couses issues in OSX

  os.remove(filename)
  os.rename(temp_filename, filename)


def patch_control(model_name, patch_dir):
  """
  Patches $SPheno/src/Control.F90
  """
  
  filename = patch_dir + "/src/Control.F90"
  temp_filename = filename + "_temp"

  # TODO: throw some errors if the file does not exist

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

  # TODO: throw some errors if the file does not exist

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

  # TODO: throw some errors if the file does not exist

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

  # TODO: throw some errors if the file does not exist

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

  # TODO: throw some errors if the file does not exist

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
  channels = {"Cha": {"ChaToChacChaCha", "ChaToChaChiChi"},
              "Chi": {"ChiToChicChaCha", "ChiToChiChiChi"},
              "Glu": {},
              "Sd": {"SdToChaGluSu", "SdToSdChacCha", "SdToSdChiChi", "SdToChiGluSd", "SdToGluGluSd"},
              "Su": {"SuToSuChiChi", "SuToChiGluSu", "SuToSdChicCha", "SuToGluGluSu", "SuToGluSdcCha", "SuToSuChacCha"},
              "Se": {"SeToSvChaChi", "SeToSeChacCha", "SeToSeChiChi"},
              "Sv": {"SvToSvChiChi", "SvToSeChicCha", "SvToSvChacCha"}}

  for particle in particles :
 
    filename = patch_dir + "/" + model_name + "/3-Body-Decays/" + particle + "_" + model_name + ".f90"
    temp_filename = filename + "_temp"

    # TODO: throw some errors if the file does not exist

    with open(filename, 'r') as f, open(temp_filename, 'w') as g :
      for line in f :
        if line.startswith("Use ThreeBodyPhaseSpace") :
          g.write(line)
          g.write("Use Model_Data_" + model_name + " ! Added by GAMBIT\n")
        elif any([line.startswith("Call " + channel) for channel in channels[particle]]) :
          g.write("If (CalcSUSY3BodyDecays) Then ! Added by GAMBIT\n")
          g.write(line)
        elif any([line.startswith("g" + channel + "(i_run,:,:,:) = g" + channel + "i") for channel in channels[particle]]) :
          g.write("End If ! Added by GAMBIT\n\n")
          g.write(line)
        else :
          g.write(line)

    os.remove(filename)
    os.rename(temp_filename, filename)


"""
FRONTEND ROUTINES
"""
