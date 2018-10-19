"""
Routines for writing CalcHEP (.mdl) files from MadGraph (.ufo) files.

The user can either compare between .ufo and .mdl files to check everything
is consistent between all files:
    
    python ufo_to_mdl <madgraph_folder> <calchep_folder>
    
Or just run on a set of .ufo files and generate .mdl files from scratch:

    python ufo_to_mdl <madgraph_folder>
    
Author: Sanjay Bloor (sanjay.bloor12@imperial.ac.uk)
If you use this script as a standalone please cite the GUM manual:
      arXiv:19XX.YYYYY

"""
    
import sys
import os
import re

class Particle():
    """
    Internal object representing a particle.
    """
    
    def __init__(self, mgname, pdg_code, name, antiname, spin, 
                colour, mass, width):
                
        self.mgname = mgname
        self.pdg_code = pdg_code
        self.name = name
        self.antiname = antiname
        self.spin = spin
        self.colour = colour
        self.mass = mass
        self.width = width
        
        if (name == antiname):
            self.selfconj = True
        else:
            self.selfconj = False
  
class Parameter():
    """
    Internal object representing a parameter.
    """
    
    def __init__(self, name, nature, type, value):
    
        self.name = name
        self.nature = nature
        self.type = type
        self.value = value

def check_mg_files(location):
    """
    Checks all MadGraph files are where they should be.
    """
    
    required = ['lorentz.py', 'parameters.py', 'particles.py', 'vertices.py']
    
    if os.path.exists(location):
    
        files = [f for f in os.listdir(location) if f.endswith('.py')]
        
        # Check that all needed files are present in the directory
        if set(required) <= set(files):
            pass
        else:
            print("ERROR: Not all required MadGraph files are here: please "
                   "check the directory specified.")
            sys.exit()

    
def import_mg_params(location):
    """
    Import MadGraph5 parameters from a given file location.
    """
    
    params = []
    
    # Read all as one big long string.
    with open(location + '/parameters.py', 'r') as f:
        lines = f.read().replace('\n','')
            
    # Find all matches to "Parameter(...)" -- ditch the texname.
    pars = re.findall(r'Parameter\((.*?) texname (.*?)\)', lines)
    
    for i in range(len(pars)):
        p = pars[i][0].split(', ')
        name = p[0].strip().split("=")[1].strip().strip("'")
        nature = p[1].strip().split("=")[1].strip().strip("'")
        type = p[2].strip().split("=")[1].strip().strip("'")
        value = p[3].strip().split("=")[1].strip().strip("'")
        
        if name == 'ZERO':
            continue

        par = Parameter(name, nature, type, value)
        params.append(par)
     
    return params
          
def import_mg_parts(location):
    """
    Import MadGraph5 particles from a given file location.
    """
    
    parts = []
    conjugates = []
    
    start_num = 0
    
    # Find last instance of "import"
    with open(location + '/particles.py', 'r') as f:
        for num, line in enumerate(f, 1):
            if re.match(r'import (.*?)', line): 
                start_num = num
    
    # Now read in properly
    with open(location + '/particles.py', 'r') as f:
        lines = ''.join(f.readlines()[start_num:]).replace('\n','')
                  
    # Find all matches to "X = Particle(...)"
    particles = re.findall(r'(.*?) = Particle\((.*?)\)', lines)
    
    for i in range(len(particles)): 
        mgname = particles[i][0]
                
        if "anti()" in mgname:
            anti, mgname = mgname.split('()')
            anti+= '()'
            conjugates.append(anti)            
                        
        p = particles[i][1].split(', ')
        pdg_code = p[0].strip().split("=")[1].strip().strip("'")
        name = p[1].strip().split("=")[1].strip().strip("'")
        antiname = p[2].strip().split("=")[1].strip().strip("'")
        spin = p[3].strip().split("=")[1].strip().strip("'") #2S+1
        color = p[4].strip().split("=")[1].strip().strip("'") # 1/3/8.
        mass = p[5].strip().split("=")[1].strip().strip("'") 
        width = p[6].strip().split("=")[1].strip().strip("'") 
        
        par = Particle(mgname, pdg_code, name, antiname, spin, color, mass, width)
        parts.append(par)
        
    # Go through conjugates and add them as mgantiname     
        
    return parts

def import_mg_verts(location):
    """
    Import MadGraph5 vertices from a given file location.
    """
    

def check_ch_files(location):
    """
    Check all CalcHEP files are where they should be.
    """
    
    required = ['func1.mdl', 'vars1.mdl', 'prtcls1.mdl', 'lgrng1.mdl']
    
    if os.path.exists(location):
    
        files = [f for f in os.listdir(location) if f.endswith('.mdl')]
        
        # Check that all needed files are present in the directory
        if set(required) <= set(files):
            pass
        else:
            print("ERROR: Not all required CalcHEP files are here: please "
                   "check the directory specified.")
            sys.exit()
   
def import_ch_params(location):
    """
    Import CalcHEP parameters from a given file location.
    """
    
    
def import_ch_parts(location):
    """
    Import CalcHEP particles from a given file location.
    """
    
    
def import_ch_verts(location):
    """
    Import CalcHEP vertices from a given file location.
    """
    
    
def convert(mg_location):
    """
    Create CalcHEP files from a given set of MadGraph files.
    """
    
    check_mg_files(mg_location)
    params = import_mg_params(mg_location)
    parts = import_mg_parts(mg_location)
    verts = import_mg_verts(mg_location)
    
    
def compare(mg_location, ch_location):
    """
    Compare a list of vertices from MadGraph with those from CalcHEP.
    By this, we simply check that each vertex is present in each set of 
    model files.
    """
    
    check_mg_files(mg_location)
    check_ch_files(ch_location)
    
    params = import_mg_params(mg_location)
    parts = import_mg_parts(mg_location)
    verts = import_mg_verts(mg_location)
    
def usage():
    """
    Print usage message.
    """
    
    msg = (
        "Usage:\n\n"
        "To compare between a set of .ufo and .mdl files, do:\n\n"
        "\tpython ufo_to_mdl <madgraph_folder> <calchep_folder>\n\n"
        "To run on a set of .ufo files and generate .mdl files, do:\n\n"
        "\tpython ufo_to_mdl <madgraph_folder>\n"
    )
    print(msg)


if __name__ == '__main__':

    print("Welcome to the .ufo to .mdl file converter!")
    print("If you use this, please cite the GUM manual: arXiv:19XX.YYYYY.\n\n")
    
    if len(sys.argv) == 1:
        usage()
        
    elif len(sys.argv) == 2:
        convert(sys.argv[1])
        
    elif len(sys.argv) == 3:
        compare(sys.argv[1], sys.argv[2])
