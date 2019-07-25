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

    patch_main_dir = patch_dir + "/"
    patch_src_dir = patch_dir + "/src/"
    patch_model_dir = patch_dir + "/" + model_name + "/"

    patch_spheno_makefile(model_name, patch_main_dir)
    patch_spheno_model_makefile(model_name, patch_model_dir)
    patch_spheno_src_makefile(model_name, patch_src_dir)
    #patch_control(model_name, patch_dir)
    #patch_brs(model_name, patch_dir)
    #patch_loopfunctions(model_name, patch_dir)
    #patch_spheno_model(model_name, patch_dir)
    #patch_model_data(model_name, patch_dir)

    # TODO: if gum.is_susy: ...
    #patch_3_body_decays_susy(model_name, patch_dir)

    print("SPheno files patched.")

# Model-independent patches

def patch_spheno_makefile(model_name, patch_dir):
  """
  Patches $SPheno/Makefile
  """

  filename = patch_dir + "Makefile"

  makefile_content = "# please put here your preferred F95/F2003 compiler\n"\
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

  f = open(filename, 'w')
  f.write(makefile_content)

def patch_spheno_model_makefile(model_name, patch_dir):
	"""
	Patches $SPheno/<MODEL>/Makefile
	"""

def patch_spheno_src_makefile(model_name, patch_dir):
	"""
	Patches $SPheno/src/Makefile
	"""

def patch_control(model_name, patch_dir):
	"""
	Patches $SPheno/src/Control.f90
	"""

def patch_brs(model_name, patch_dir):
	"""
	Patches $SPheno/BranchingRatios_<MODEL>.f90
	"""

def patch_loopfunctions(model_name, patch_dir):
	"""
	Patches $SPheno/AddLoopFunctions.f90
	"""

# Model-dependent patches

def patch_spheno_model(model_name, patch_dir):
	"""
	Patches $SPheno/SPheno<MODEL>.f90
	"""

	filename = "{0}/SPheno{1}.f90".format(patch_dir, model_name)
	tempfile = filename+"_temp"

	with open(filename, 'r') as f, open(tempfile, 'w') as g:
		for line in f:
			if line.startswith("Program SPheno"+model_name):
				g.write("!Program SPheno"+model_name+" ! Commented by GAMBIT")
				g.write("Module SPheno"+model_name+" ! Added by GAMBIT")
			elif line.startswith("End Program SPheno"+model_name):
				g.write("!End Program SPheno"+model_name+" ! Commented by GAMBIT")
				g.write("End Module SPheno"+model_name+" ! Added by GAMBIT")
			else:
				g.write(line)


def patch_model_data(model_name, patch_dir):
	"""
	Patches $SPheno/Model_Data_<MODEL>.f90
	"""

# SUSY-only patches

def patch_3_body_decays_susy(model_name, patch_dir):
	"""
	Patches the 3-body decays in: 
	$SPheno/3-Body-Decays/X_<MODEL>.f90
	where X is a superfield.
	"""

"""
FRONTEND ROUTINES
"""
