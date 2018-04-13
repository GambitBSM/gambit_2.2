"""
Master module containing all routines for finding, creating,
amending, and reading files.
"""

import os
import re
import numpy as np
from setup import *

mode = 'Test'
#mode = 'Go'

def full_filename(filename, module, header=False):
    """
    Formats a gambit file correctly based on the filename, the 
    module, and whether it is a header file.
    """
    
    # strip leading & trailing slashes
    filename.strip('/')
    module.strip('/')
    if not module.endswith('Bit') and module != "Models":
        module += 'Bit'
    
    path = ""
    if header == True:
        filename += ".hpp"
        path = "include/gambit/{0}".format(module)
    else:
        filename += ".cpp"
        path = "src"
    
    location = "../{0}/{1}/{2}".format(module, path, filename)
    return location

def find_file(filename, module, header=False):
    """
    Tries to find a file in a specified module.
    """
    
    location = full_filename(filename, module, header)
    
    if os.path.exists(location):
        return True
    else:
        return False
    
def find_capability(capability, module, filename=None):
    """
    Tries to find a capability in a specified rollcall file.
    """
    
    module.strip('/')
    if not module.endswith('Bit') and module != "Models":
        module += 'Bit'
    
    if not filename:
        filename = "{0}_rollcall".format(module)
    else:
        filename.strip('/')
  
    location = full_filename(filename, module, header=True)
    if find_file:
        pass
    else:
        GumError(("Cannot find capability {0} in rollcall header file {1} "
                  "as the file does not exist!!").format(capability, 
                                                         location))
  
    lookup = "#define CAPABILITY " + capability
    
    with open(location) as f:
        for num, line in enumerate(f, 1):
            if lookup in line:
                return True, num
          
    return False, 0
  
def find_function(function, capability, module, filename=None):
    """
    Tries to find a function in a specified rollcall header file.
    """
    
    module.strip('/')
    if not module.endswith('Bit') and module != "Models":
        module += 'Bit'
    
    if not filename:
        filename = "{0}_rollcall".format(module)
    else:
        filename.strip('/')
    
    location = full_filename(filename, module, header=True)
  
    # First check the capability exists...
    exists, num = find_capability(capability, module, filename)
    
    lookup = "#define FUNCTION " + function
    terminate = "#undef CAPABILITY"
    
    with open(location) as f:
        for i in xrange(num):
            f.next()
        for no, line in enumerate(f, 1+num):
            if lookup in line:
                return True, no
            if terminate in line:
                return False, 0
    
    return False, 0
  
def write_file(filename, module, contents, header=False):
    """
    Writes a file in a specified location.
    """
    
    location = full_filename(filename, module, header)
      
    if find_file(filename, module, header) == True:
        raise GumError(("ERROR: Tried to write file " + location + 
                        ", but it already exists."))
                      
    if mode != 'Test':
        # Create new file 
        open(location, 'w').write(contents)
  
    print("File {} successfully created.".format(location))
  
def delete_file(filename, module, header=False):
    """
    Deletes a file in a specified location.
    """
    
    location = full_filename(filename, module, header)
    
    if find_file(filename, module, header) == True:
        os.remove(location)
        print("File {} successfully removed.".format(location))

def amend_file(filename, module, contents, line_number, header=False):
    """
    Amends a file in a specified location with 'contents', starting 
    from a given line number.
    """
    
    location = full_filename(filename, module, header)
    
    if find_file(filename, module, header) == False:
        raise GumError(("ERROR: Tried to amend file " + location + 
                        ", but it does not exist."))
    
    if mode != 'Test':
        temp_location = location + "_temp"
        lines = open(location, 'r').readlines()
      
        with open(temp_location, 'w') as f:
            for i in xrange(line_number):
                f.write(lines[i])
            f.write(contents)
            for i in xrange(len(lines)-line_number):
                f.write(lines[i+line_number])
      
      os.remove(location)
      os.rename(temp_location, location)

  print("File {} successfully amended.".format(location))
  
def write_capability(capability, functions):
    """
    Writes a new capability, wrapped around (a) function(s), 
    for a rollcall header file.
    """
      
    towrite = (
            "  #define CAPBILITY {0}\n"
            "  START_CAPABILITY\n"
            "  \n"
            "{1}"
            "  \n"
            "  #undef CAPABILITY\n"
    ).format(capability, '\n'.join(functions))
    
    print("Capability {} added.".format(capability))
    
    return towrite

def write_function(function, returntype, dependencies=None, 
                   allowed_models=None, backend_reqs=None):
    """
    Writes a function for a rollcall header file.
    """
    
    extras = ""
    
    if dependencies:
        for i in np.arange(len(dependencies)):
            extras += "DEPENDENCY({0}, {1})\n".format(dependencies[i][0], 
                                                      dependencies[i][1])
    if allowed_models:
        extras += "ALLOW_MODELS({0})\n".format(', '.join(allowed_models))
      
    if backend_reqs:
        for i in np.arange(len(backend_reqs)):
            bereq = backend_reqs[i]
            if not isinstance(bereq, BackendReq):
                raise GumError(("Backend requirement at position " + i + 
                                "not passed as instance of class BackendReq."))
                          
            add = ", ".join([bereq.capability, 
                             "({0})".format(', '.join(bereq.tags)), 
                             bereq.cpptype])
        
            if not bereq.var:
                add += ", ({0})".format(', '.join(bereq.args))     
                      
            extras += "BACKEND_REQ({0})\n".format(add)
      
    towrite = (
            "\n"
            "#define FUNCTION {0}\n"
            "START_FUNCTION({1})\n"
            "{2}"
            "#undef FUNCTION\n"
    ).format(function, returntype, extras)
    
    return dumb_indent(4, towrite)
  
