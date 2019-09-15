"""
A dummy example of a python code providing 
LHC cross-sections to GAMBIT  
"""
from __future__ import print_function
import xsec 
import numpy as np 


#
# This code will run when the library is loaded
#
prefix = "xsec_example_gambit:"
print(prefix)
print(prefix, "Starting up...")
print(prefix, "Have loaded xsec version", xsec.__version__, "from", xsec.__path__[0])
print(prefix, "Have loaded numpy version", np.__version__, "from", np.__path__[0])


#
# Set global parameters
#
def set_parameters(params):
    prefix = "xsec_example_gambit: set_parameters:"

    print(prefix)
    print(prefix, "Starting function...")
    print(prefix, "I got a 'params' dictionary with these key-value pairs:")
    for k,v in params.items():
        print(prefix, "-", k, ":", v)


#
# Set global library flags
#
def set_flags(flags):
    prefix = "xsec_example_gambit: set_flags:"

    print(prefix)
    print(prefix, "Starting function...")
    print(prefix, "I got a 'flags' dictionary with these key-value pairs:")
    for k,v in flags.items():
        print(prefix, "-", k, ":", v)



#
# Return the cross-section for a given process, identified by a PID pair
#
def xsec_fb(proc, proc_params, proc_flags):
    prefix = "xsec_example_gambit: xsec_fb:"

    print(prefix)
    print(prefix, "Starting function...")
    print(prefix, "I got this 'proc' array:", proc)
    print(prefix)
    print(prefix, "I got a 'proc_params' dictionary with these key-value pairs:")
    for k,v in proc_params.items():
        print(prefix, "-", k, ":", v)
    print(prefix)
    print(prefix, "I got a 'proc_flags' dictionary with these key-value pairs:")
    for k,v in proc_flags.items():
        print(prefix, "-", k, ":", v)

    return 123.45


#
# Return the asymmetric cross-section error for a given process, identified by a PID pair
#
def xsec_err_fb(proc, proc_params, proc_flags):
    prefix = "xsec_example_gambit: xsec_err_fb:"

    print(prefix)
    print(prefix, "Starting function...")
    print(prefix, "I got this 'proc' array:", proc)
    print(prefix)
    print(prefix, "I got a 'proc_params' dictionary with these key-value pairs:")
    for k,v in proc_params.items():
        print(prefix, "-", k, ":", v)
    print(prefix)
    print(prefix, "I got a 'proc_flags' dictionary with these key-value pairs:")
    for k,v in proc_flags.items():
        print(prefix, "-", k, ":", v)

    return np.array([1.1, 1.2])




# # ================ TESTING ================

# params_in = {'par1': 1.0, 'par2': 2.0}
# set_parameters(params_in)

# flags_in = {'flag1': True, 'flag2': False}
# set_flags(flags_in)

# xsec_fb([1000021,1000022], {'LO_xsec_fb': 78.9}, {})

# xsec_err_fb([1000021,1000022], {}, {'alphas_error': True})




