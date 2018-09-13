from src import *

import sys
path = os.path.join(os.getcwd(),'contrib','MadGraph')
sys.path.append(path)

import madgraph.interface.master_interface as mi

launch = mi.MasterCmd(mgme_dir = path)
launch.exec_cmd("import model SingletDM_test")
launch.exec_cmd("generate process p p > ~s ~s")
launch.exec_cmd("output pythia8")
