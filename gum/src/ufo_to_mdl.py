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

"""
MadGraph objects and routines
"""

class MGParticle():
    """
    Internal object representing a particle
    within MadGraph.
    """
    
    def __init__(self, mgname, mgantiname, pdg_code, name, antiname, spin, 
                 color, mass, width):
                
        self.mgname = mgname
        self.mgantiname = mgantiname
        self.pdg_code = pdg_code
        self.name = name
        self.antiname = antiname
        self.spin = spin
        self.color = color
        self.mass = mass
        self.width = width
        
        if (name == antiname):
            self.selfconj = True
        else:
            self.selfconj = False
  
class MGParameter():
    """
    Internal object representing a parameter
    within MadGraph.
    """
    
    def __init__(self, name, nature, type, value):
    
        self.name = name
        self.nature = nature
        self.type = type
        self.value = value

class MGCoupling():
    """
    Internal object representing a coupling structure
    within MadGraph.
    """

    def __init__(self, mgname, value, order):

        self.mgname = mgname
        self.value = value
        self.order = order

class MGVertex():
    """
    Internal object representing a vertex
    within MadGraph.
    """

    def __init__(self, mgname, particles, color, lorentz, couplings):

        self.mgname = mgname
        self.particles = particles
        self.color = color
        self.lorentz = lorentz
        self.couplings = couplings
        self.numparticles = len(particles)


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

        par = MGParameter(name, nature, type, value)
        params.append(par)
     
    return params
          
def import_mg_parts(location):
    """
    Import MadGraph5 particles from a given file location.
    """
    
    parts = []
    
    start_num = 0
    
    # Find last instance of "import"
    with open(location + '/particles.py', 'r') as f:
        for num, line in enumerate(f, 1):
            if re.match(r'import (.*?)', line): 
                start_num = num
    
    # Now read in properly
    with open(location + '/particles.py', 'r') as f:
        lines = ''.join(f.readlines()[start_num:]).replace('\n','')

    # Explicitly add a newline after each closing parenthesis.
    # This makes the regular expressions work consistently.
    lines = lines.replace(')', ')\n')

    # Find all matches to "X = Particle(...)"
    particles = re.findall(r'(.*?) = Particle\((.*?)\)\n', lines)

    # Find all matches to "Y = Y'.anti()"
    antis = re.findall(r'(.*?) = (.*?).anti\(\)\n', lines)
    
    for i in range(len(particles)): 

        mgname = particles[i][0]

        p = particles[i][1].split(', ')

        pdg_code = p[0].strip().split("=")[1].strip().strip("'")
        name = p[1].strip().split("=")[1].strip().strip("'")
        antiname = p[2].strip().split("=")[1].strip().strip("'")
        spin = p[3].strip().split("=")[1].strip().strip("'")     # 2S+1
        color = p[4].strip().split("=")[1].strip().strip("'")    # 1/3/8.
        mass = p[5].strip().split("=")[1].strip().strip("'") 
        width = p[6].strip().split("=")[1].strip().strip("'") 

        # Go through conjugates and add them as mgantiname -- if distinct.
        mgantiname = mgname
        for j in range(len(antis)):
            if (antis[j][1] == mgname):
                mgantiname = antis[j][0]
        
        par = MGParticle(mgname, mgantiname, pdg_code, name, antiname, 
                       spin, color, mass, width)
        parts.append(par)

    return parts

def import_mg_couplings(location):
    """
    Import MadGraph5 couplings from a given file location.
    """

    couplings = []
    
    start_num = 0
    
    # Find last instance of "from ___ import ___"
    with open(location + '/couplings.py', 'r') as f:
        for num, line in enumerate(f, 1):
            if re.match(r'from (.*?)', line): 
                start_num = num
    
    # Now read in properly
    with open(location + '/couplings.py', 'r') as f:
        lines = ''.join(f.readlines()[start_num:]).replace('\n','')

    # Explicitly add a newline after each closing parenthesis.
    # This makes the regular expressions work consistently.
    lines = lines.replace('})', '})\n')

    # Find all matches to "X = Vertex(...)"
    couplings = re.findall(r'(.*?) = Coupling\((.*?)\)\n', lines)

    for i in range(len(couplings)):

        mgname = couplings[i][0]
        
        # Match 'value = VALUE'
        value = re.findall(r'value = \'(.*?)\'', str(couplings[i]))[0]
        order = re.findall(r'order = {(.*?)}', str(couplings[i]))[0]

        coupling = MGCoupling(mgname, value, order)
        couplings.append(coupling)

    return couplings

def import_mg_verts(location):
    """
    Import MadGraph5 vertices from a given file location.
    """

    verts = []
    
    start_num = 0
    
    # Find last instance of "import"
    with open(location + '/vertices.py', 'r') as f:
        for num, line in enumerate(f, 1):
            if re.match(r'import (.*?)', line): 
                start_num = num
    
    # Now read in properly
    with open(location + '/vertices.py', 'r') as f:
        lines = ''.join(f.readlines()[start_num:]).replace('\n','')

    # Explicitly add a newline after each closing parenthesis.
    # This makes the regular expressions work consistently.
    lines = lines.replace('})', '})\n')

    # Find all matches to "X = Vertex(...)"
    vertices = re.findall(r'(.*?) = Vertex\((.*?)\)\n', lines)

    for i in range(len(vertices)):

        v = vertices[i][1].split(', ')
        mgname = v[0].strip().split("=")[1].strip().strip("'")

        # Match 'particles = [ X, Y, Z ]'
        parts = re.findall(r'particles = \[(.*?)\]', str(vertices[i]))
        # This returns a list with 1 entry; we want to split this up.
        parts = parts[0].split(',')
        # Remove those pesky whitespaces
        particles = [x.strip() for x in parts]

        # 'Color' entry
        color = re.findall(r'color = \[(.*?)\]', str(vertices[i]))
        color = color[0].strip().strip("'")

        # 'Lorentz' entry
        lorentz = re.findall(r'lorentz = \[(.*?)\]', str(vertices[i]))
        # This returns a list with 1 entry; we want to split this up.
        lorentz = lorentz[0].split(',')
        # Remove whitespace and quote marks
        lorentz = [x.strip() for x in lorentz]

        couplings = re.findall(r'couplings = {(.*?)}', str(vertices[i]))
        
        vertex = MGVertex(mgname, particles, color, lorentz, couplings)
        verts.append(vertex)


    return verts

"""
CalcHEP objecs and routines
"""

class CHParticle():
    """
    Internal object representing a particle
    within CalcHEP.
    """
    
    def __init__(self, chname, pdg_code, name, antiname, spinx2, 
                 color, mass, width, aux):
                
        self.chname = chname
        self.pdg_code = pdg_code
        self.name = name
        self.antiname = antiname
        self.spinx2 = spinx2
        self.color = color
        self.mass = mass
        self.width = width
        self.aux = aux
        
        if (name == antiname):
            self.selfconj = True
        else:
            self.selfconj = False
  
class CHParameter():
    """
    Internal object representing a parameter
    within CalcHEP.
    """
    
    def __init__(self, name, nature, value):
    
        self.name = name
        self.nature = nature
        self.value = value

class CHCoupling():
    """
    Internal object representing a coupling structure
    within CalcHEP.
    """

    def __init__(self, chname, value, order):

        self.chname = chname
        self.value = value
        self.order = order

class CHVertex():
    """
    Internal object representing a vertex
    within CalcHEP.
    """

    def __init__(self, particles, lorentz, couplings):

        self.particles = particles
        self.lorentz = lorentz
        self.couplings = couplings
        self.numparticles = len(particles)


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
    
    params = []
    start_num = 0

    # Firstly, deal with internal parameters

    # Trim the stuff of the top
    with open(location + '/func1.mdl') as f:
        for num, line in enumerate(f, 1):
            if re.search(r'Expression', line):
                start_num = num
    
    # Now read in line by line (since all parameter info is stored on 1 line)
    with open(location + '/func1.mdl', 'r') as f:
        lines = f.readlines()[start_num:]

    for line in lines:

        # Get all information in a nice list, remove the whitespace
        p = [i.strip() for i in line.split('|')]

        name = p[0]

        # FeynRules output comments out its SPheno interface (!) if you don't request it
        # - don't want to parse this.
        if name.startswith('%'):
            continue

        # All of these are internal.
        nature = "Internal"

        # Remove any comments from the value
        value = p[1].split('%')[0].strip()

        parameter = CHParameter(name, nature, value)
        params.append(parameter)

    # Now external parameters

    with open(location + '/vars1.mdl') as f:
        for num, line in enumerate(f, 1):
            if re.search(r'Value', line):
                start_num = num

    with open(location + '/func1.mdl', 'r') as f:
        lines = f.readlines()[start_num:]

    for line in lines:

        p = [i.strip() for i in line.split('|')]

        name = p[0]
        value = p[1]
        nature = "External"

        parameter = CHParameter(name, nature, value)
        params.append(parameter)

    return params
    
def import_ch_parts(location):
    """
    Import CalcHEP particles from a given file location.
    """

    parts = []
    
    start_num = 0
    
    # Trim the stuff of the top
    with open(location + '/prtcls1.mdl') as f:
        for num, line in enumerate(f, 1):
            if re.search(r'Full  name', line):
                start_num = num
    
    # Now read in line by line (since all particle info is stored on 1 line)
    with open(location + '/prtcls1.mdl', 'r') as f:
        lines = f.readlines()[start_num:]

    for line in lines:

        # Get all information in a nice list, remove the whitespace
        p = [i.strip() for i in line.split('|')]
        
        chname = p[0]
        name = p[1]
        antiname = p[2]
        pdg_code = p[3]
        spinx2 = p[4]
        mass = p[5]
        width = p[6]
        color = p[7]
        aux = p[8]

        particle = CHParticle(chname, pdg_code, name, antiname, spinx2, color, mass, width, aux)
        parts.append(particle)

    return parts

    
def import_ch_verts(location):
    """
    Import CalcHEP vertices from a given file location.
    """

    verts = []
    start_num = 0

    # Trim the stuff of the top
    with open(location + '/lgrng1.mdl') as f:
        for num, line in enumerate(f, 1):
            if re.search(r'Factor', line):     
                start_num = num

    # Now read in line by line (since all vertex info is stored on 1 line)
    with open(location + '/lgrng1.mdl', 'r') as f:
        lines = f.readlines()[start_num:]

    for line in lines:

        # Get all information in a nice list, remove the whitespace
        v = [i.strip() for i in line.split('|')]

        # Always 3 particles in a vertex...
        particles = [v[0], v[1], v[2]]
        # Sometimes a fourth - if not empty, add it
        if v[3] != '':
            particles.append(v[3])
        
        coupling = v[4]
        lorentz = v[5]

        vertex = CHVertex(particles, lorentz, coupling)
        verts.append(vertex)

    return verts

def convert(mg_location):
    """
    Create CalcHEP files from a given set of MadGraph files.
    """
    
    # Check the MadGraph files are where the user says
    check_mg_files(mg_location)

    # Read in the MadGraph files
    params = import_mg_params(mg_location)
    parts = import_mg_parts(mg_location)
    couplings = import_mg_couplings(mg_location)
    verts = import_mg_verts(mg_location)
    
    
def compare(mg_location, ch_location):
    """
    Compare a list of vertices from MadGraph with those from CalcHEP.
    By this, we simply check that each vertex is present in each set of 
    model files.
    """
    
    # Check both sets of files are where they should be
    check_mg_files(mg_location)
    check_ch_files(ch_location)
    
    # Read in the MadGraph files
    mg_params = import_mg_params(mg_location)
    mg_parts = import_mg_parts(mg_location)
    mg_couplings = import_mg_couplings(mg_location)
    mg_verts = import_mg_verts(mg_location)

    # Read in the CalcHEP files
    ch_params = import_ch_params(ch_location)
    ch_parts = import_ch_parts(ch_location)
    ch_verts = import_ch_verts(ch_location)
    
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
