from src import *

import sys

# Add MadGraph to path
path = os.path.join(os.getcwd(),'contrib','MadGraph')
sys.path.append(path)

# Go to the MadGraph directory so we don't clutter up $GUM
os.chdir(path)

import madgraph.interface.master_interface as mi

print("Launching MadGraph...")

# Launch interactive interface to MadGraph
launch = mi.MasterCmd(mgme_dir = path)

# Import the model
launch.exec_cmd("import model SingletDM_test")

# Import *all* cross-sections we want to generate
launch.exec_cmd("generate p p > ~s ~s")

# Next ones have the format
# launch.exec_cmd("add process p p > u u~")

# Create output for Pythia
launch.exec_cmd("output pythia8")

print("Pythia output created.")
