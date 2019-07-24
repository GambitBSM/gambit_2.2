"""
Contains all routines for SPheno output from SARAH.
"""

from files import *
from distutils.dir_util import copy_tree
from collections import defaultdict
import re

# An empty defaultdict in scope, so we don't save to the 
# mug file, but can still use the writing routines in 
# the 'files' module.
d = defaultdict(list)

def copy_spheno_files(model_name, output_dir, spheno_path):
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

		# Now copy from SARAH to the Output directory
		copy_tree(spheno_path, newdir)

	print("SPheno files moved to output directory, and a copy made.")

"""
PATCHING
"""

def patch_spheno(model_name, patch_dir):
	"""
	Applies all patches to SPheno in the GUM 
	Outputs/... directory.
	"""

	# SARAH removes hyphens for .f90 files...
	clean_model_name = model_name.replace('-','') 

	patch_spheno_makefile(clean_model_name, patch_dir)
	patch_spheno_model_makefile(clean_model_name, patch_dir)
	patch_spheno_src_makefile(clean_model_name, patch_dir)
	patch_control(clean_model_name, patch_dir)
	patch_brs(clean_model_name, patch_dir)
	patch_loopfunctions(clean_model_name, patch_dir)
	patch_spheno_model(clean_model_name, patch_dir)
	patch_model_data(clean_model_name, patch_dir)

	# TODO: if gum.is_susy: ...
	patch_3_body_decays_susy(clean_model_name, patch_dir)

	print("SPheno files patched.")

# Model-independent patches

def patch_spheno_makefile(clean_model_name, patch_dir):
	"""
	Patches $SPheno/Makefile
	"""

def patch_spheno_model_makefile(clean_model_name, patch_dir):
	"""
	Patches $SPheno/<MODEL>/Makefile
	"""

def patch_spheno_src_makefile(clean_model_name, patch_dir):
	"""
	Patches $SPheno/src/Makefile
	"""

def patch_control(clean_model_name, patch_dir):
	"""
	Patches $SPheno/src/Control.f90
	"""

def patch_brs(clean_model_name, patch_dir):
	"""
	Patches $SPheno/BranchingRatios_<MODEL>.f90
	"""

def patch_loopfunctions(clean_model_name, patch_dir):
	"""
	Patches $SPheno/AddLoopFunctions.f90
	"""

# Model-dependent patches

def patch_spheno_model(clean_model_name, patch_dir):
	"""
	Patches $SPheno/SPheno<MODEL>.f90
	"""

	filename = "{0}/SPheno{1}.f90".format(patch_dir, clean_model_name)
	tempfile = filename+"_temp"

	with open(filename, 'r') as f, open(tempfile, 'w') as g:
		for line in f:
			if line.startswith("Program SPheno"+clean_model_name):
				g.write("!Program SPheno"+clean_model_name+" ! Commented by GAMBIT")
				g.write("Module SPheno"+clean_model_name+" ! Added by GAMBIT")
			elif line.startswith("End Program SPheno"+clean_model_name):
				g.write("!End Program SPheno"+clean_model_name+" ! Commented by GAMBIT")
				g.write("End Module SPheno"+clean_model_name+" ! Added by GAMBIT")
			else:
				g.write(line)


def patch_model_data(clean_model_name, patch_dir):
	"""
	Patches $SPheno/Model_Data_<MODEL>.f90
	"""

# SUSY-only patches

def patch_3_body_decays_susy(clean_model_name, patch_dir):
	"""
	Patches the 3-body decays in: 
	$SPheno/3-Body-Decays/X_<MODEL>.f90
	where X is a superfield.
	"""

"""
FRONTEND ROUTINES
"""