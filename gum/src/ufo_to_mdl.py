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
import shutil

counter = 0

"""
MadGraph objects and routines
"""

class MGParticle():
    """
    Internal object representing a particle
    within MadGraph.
    """
    
    def __init__(self, mgname, mgantiname, pdg_code, name, antiname, spin, 
                 color, mass, width, goldstone, ghostnum):
                
        self.mgname = mgname
        self.mgantiname = mgantiname
        self.pdg_code = pdg_code
        self.name = name
        self.antiname = antiname
        self.spin = spin
        self.color = color
        self.mass = mass
        self.width = width
        self.goldstone = goldstone
        self.ghostnum = ghostnum
        
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

class MGLorentz():
    """
    Internal object representing a lorentz structure
    within MadGraph.
    """

    def __init__(self, mgname, spins, structure):

        self.mgname = mgname
        self.spins = spins   #2s+1
        self.structure = structure

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

        # Only care about external parameters (for now).
        if nature == 'external': params.append(par)
    
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

        p = particles[i][1]

        pdg_code = re.findall(r'pdg_code = (.*?),', p)[0].strip("'")
        name = re.findall(r'name = (.*?),', p)[0].strip("'")
        antiname = re.findall(r'antiname = (.*?),', p)[0].strip("'")
        spin = re.findall(r'spin = (.*?),', p)[0].strip("'")     # 2S+1
        color = re.findall(r'color = (.*?),', p)[0].strip("'")    # 1/3/8.
        mass = re.findall(r'mass = (.*?),', p)[0].strip("'") 
        width = re.findall(r'width = (.*?),', p)[0].strip("'")

        # Also save the ghost number, they're represented differently.
        ghostnum = re.findall(r'GhostNumber = (.*?),', p)[0]

        # Same with Goldstones
        goldstone = re.findall(r'goldstone = (.*?),', p)
        if len(goldstone) > 0: goldstone = goldstone[0]
        else: goldstone = False

        # Go through conjugates and add them as mgantiname -- if distinct.
        mgantiname = mgname
        for j in range(len(antis)):
            if (antis[j][1] == mgname):
                mgantiname = antis[j][0]
        
        par = MGParticle(mgname, mgantiname, pdg_code, name, antiname, 
                       spin, color, mass, width, goldstone, ghostnum)
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

def import_mg_lorentz(location):
    """
    Import Lorentz structures from MadGraph5 from a given file location.
    """

    lor = []

    start_num = 0

    # Find last instance of "pass"
    with open(location + '/lorentz.py', 'r') as f:
        for num, line in enumerate(f, 1):
            if re.match(r'from (.*?)', line): 
                start_num = num

    # Now read in properly
    with open(location + '/lorentz.py', 'r') as f:
        lines = ''.join(f.readlines()[start_num:]).replace('\n','')

    # Explicitly add a newline after each closing parenthesis.
    # This makes the regular expressions work consistently.
    lines = lines.replace('\')', '\')\n')

    # Find all matches to "X = Particle(...)"
    lorentz = re.findall(r'(.*?) = Lorentz\((.*?)\)\n', lines)

    for i in range(len(lorentz)):

        l = lorentz[i][1]

        # Match 'name = '...' '
        mgname = re.findall(r'name = \'(.*?)\'', str(l))[0]

        # Match 'spins = [ X, Y, Z ]'
        spins = re.findall(r'spins = \[(.*?)\]', str(l))
        # This returns a list with 1 entry; we want to split this up.
        spins = spins[0].split(',')
        # Remove those pesky whitespaces
        spins = [x.strip() for x in spins]

        # Match 'structure = '...' '
        structure = re.findall(r'structure = \'(.*?)\'', str(l))[0]
        
        lorentztry = MGLorentz(mgname, spins, structure)
        lor.append(lorentztry)

    return lor

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
        particles = [x.strip()[2:] for x in parts]

        # 'Color' entry
        color = re.findall(r'color = \[(.*?)\]', str(vertices[i]))[0]
        # We have a list of color indices now, all between 
        # inverted commas so match these
        color = re.findall(r'\'(.*?)\'', color)
        # Create a list object
        color = [x.strip() for x in color]

        # 'Lorentz' entry
        lorentz = re.findall(r'lorentz = \[(.*?)\]', str(vertices[i]))
        # This returns a list with 1 entry; we want to split this up.
        lorentz = lorentz[0].split(',')
        # Remove whitespace and quote marks
        lorentz = [x.strip() for x in lorentz]

        couplings = re.findall(r'couplings = {(.*?)}', str(vertices[i]))

        # Commas appear in the coupling structure -- use the parenthesis too
        couplings = couplings[0].split(",(")

        # If there were more than 1, add the parentheses back
        if len(couplings) > 1:
            for i in range(len(couplings)-1):
                couplings[i+1] = '(' + couplings[i+1]

        # Each coupling should have a Lorentz structure with it...
        if len(lorentz) != len(couplings):
            print("Length of Lorentz entry not the same as couplings entry!")
            print("The vertex that caused the problem is:")
            print(vertices[i][1])
            sys.exit()
        
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


    # Currently only deal with external parameters

    # # Trim the stuff off the top
    # with open(location + '/func1.mdl') as f:
    #     for num, line in enumerate(f, 1):
    #         if re.search(r'Expression', line):
    #             start_num = num
    
    # # Now read in line by line (since all parameter info is stored on 1 line)
    # with open(location + '/func1.mdl', 'r') as f:
    #     lines = f.readlines()[start_num:]

    # for line in lines:

    #     # Get all information in a nice list, remove the whitespace
    #     p = [i.strip() for i in line.split('|')]

    #     name = p[0]

    #     # FeynRules output comments out its SPheno interface (!) if you don't request it
    #     # - don't want to parse this.
    #     if name.startswith('%'):
    #         continue

    #     # All of these are internal.
    #     nature = "internal"

    #     # Remove any comments from the value
    #     value = p[1].split('%')[0].strip()

    #     parameter = CHParameter(name, nature, value)
    #     params.append(parameter)

    # Now external parameters

    with open(location + '/vars1.mdl') as f:
        for num, line in enumerate(f, 1):
            if re.search(r'Value', line):
                start_num = num

    with open(location + '/vars1.mdl', 'r') as f:
        lines = f.readlines()[start_num:]

    for line in lines:

        p = [i.strip() for i in line.split('|')]

        name = p[0]
        value = p[1]

        # FeynRules output comments out its SPheno interface (!) if you don't request it
        # - don't want to parse this.
        if name.startswith('%'):
            continue

        # Also ignore if it's exp or pi
        if name == 'Pi' or name == 'E':
            continue

        nature = "external"

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

"""
INTERFACING ROUTINES
"""

def get_mg_ch_param_dict():
    """
    Dictionary of common (Standard Model) parameter names
    that systematically differ between MadGraph and CalcHEP.

    The format is (key, value) = (MG, CH)
    """

    param_dict = {}

    
    param_dict["G"]         = "GG"     # Strong coupling 
    param_dict["ee"]        = "EE"     # EM constant

    return param_dict

def get_mg_ch_lorentz_dict():
    """
    Returns a dictionary to translate between Lorentz structures in MadGraph
    and CalcHEP.
    """

    ld = {}
    return ld

def get_mg_ch_parts_dict(mg_parts, ch_parts):
    """
    Dictionary of names of particles in MadGraph and CalcHEP, 
    for writing routines.

    The format is (key, value) = (MG, CH)
    """

    # Get a list of all PDG codes according to MG (ignoring ghosts & Goldstones)
    mg_pdgs = [x.pdg_code for x in mg_parts if int(x.ghostnum) == 0 and x.goldstone == False]

    ch_pdgs = [x.pdg_code for x in ch_parts]

    # Check there's no particles missing from either
    if (set(mg_pdgs) ^ set(ch_pdgs)):
        print("Not all particles in MadGraph are present in CalcHEP!")
        print("The following particles by PDG code are in one but not the other:")
        print(set(mg_pdgs) ^ set(ch_pdgs))
        print("Please check your files.")
        sys.exit()

    # We're all good now. Now go through the list of MadGraph particles, and get their names
    mg_pairs = {}
    for i in range(len(mg_parts)):
        mg_pairs[mg_parts[i].pdg_code] = mg_parts[i].name

    keys = []
    values = []

    for i in range(len(ch_parts)):
        keys.append( mg_pairs.get(ch_parts[i].pdg_code) )
        values.append( ch_parts[i].name)

    # Now construct a dictionary of {'MadGraph name': 'CalcHEP name'}
    part_dict = dict(zip(keys, values))

    return part_dict

def run_through_parameter_dict(parameters, param_dict):
    """
    Goes through all unshared parameters to check which ones 
    are not duped, just have a different name.
    """

    chparams = []
    mgparams = []

    for i in range(len(parameters)):

        p = parameters[i]

        # If it's a MadGraph vertex, look to see if the key is in the dictionary
        if isinstance(p, MGParameter):
            if p.name in param_dict.keys():
                continue
            else:
                mgparams.append(p)

        # If it's a CalcHEP one, look to see if the value is in the dictionary
        elif isinstance(p, CHParameter):
            if p.name in param_dict.values():
                continue
            else:
                chparams.append(p)

        else:
            print("Parameter passed over not an instance of either the CalcHEP or MadGraph")
            print("parameter class!")
            sys.exit()

    return mgparams, chparams

def check_for_optimisation(mgparams, chparams):
    """
    Checks the definitions of the parameters...
    """

def compare_vertices(mg_verts, mg_parts, mg_params,
                     ch_verts, ch_parts, ch_params,
                     ch_location, 
                     param_dict, lorentz_dict, part_dict,
                     ):
    """
    Compares the vertices from a MadGraph model with those
    in a CalcHEP model.

    todo: compare the vertex factors + lorentz structures
    """

    # Display all vertices by PDG code for ease
    mg_pdg_dict = {}
    ch_pdg_dict = {}

    # Have a separate version for ghosts + goldstones as they're represented somewhat
    # differently between CH and MG
    mg_ghosts = {}
    ch_ghosts = {}

    for i in range(len(mg_parts)):
        p = mg_parts[i]

        if (int(p.ghostnum) == 1) or (not p.goldstone == False):
            mg_ghosts[p.mgname] = int(p.pdg_code)
            if not p.selfconj: mg_ghosts[p.mgantiname] = -1*int(p.pdg_code)
            
        else:
            mg_pdg_dict[p.mgname] = int(p.pdg_code)
            if not p.selfconj: mg_pdg_dict[p.mgantiname] = -1*int(p.pdg_code)

    for i in range(len(ch_parts)):
        p = ch_parts[i]

        ch_pdg_dict[p.name] = int(p.pdg_code)
        if not p.selfconj: ch_pdg_dict[p.antiname] = -1*int(p.pdg_code)

    # Save each individual vertex as a set of PDG codes from MadGraph
    mg_vertices_by_pdg = []
    mg_vertices_with_ghosts = []

    for i in range(len(mg_verts)):
        v = mg_verts[i]

        parts = v.particles

        # If there's ghosts/goldstones, put them in a different list
        # and deal with them separately.
        if any( part in mg_ghosts for part in parts):
            mg_vertices_with_ghosts.append(parts)
            continue

        vpdg = [ mg_pdg_dict[i] for i in parts ]

        mg_vertices_by_pdg.append(vpdg)

    # Same for CalcHEP
    ch_vertices_by_pdg = []
    ch_vertices_with_ghosts = []

    for i in range(len(ch_verts)):
        v = ch_verts[i]

        parts = v.particles

        # Find ghosts and goldstones
        # TODO: tensor structures...
        r = re.compile(r'(.*?)\.(.*?)')

        hasghosts = False

        for part in parts:
            if r.match(part):
                hasghosts = True
        if hasghosts: continue

        vpdg = [ ch_pdg_dict[i] for i in v.particles ]

        ch_vertices_by_pdg.append(vpdg)

    
    # Reorder the PDG codes in ascending order so we can compare
    ch_vertices_by_pdg = [sorted(i) for i in ch_vertices_by_pdg]
    mg_vertices_by_pdg = [sorted(i) for i in mg_vertices_by_pdg]

    # Ditch any MG vertices with more than 4 particles, as there's no way CalcHEP 
    # can possibly have these.
    mg_set = set(tuple(i) for i in mg_vertices_by_pdg if len(i) < 5)
    # Now compare them
    ch_set = set(tuple(i) for i in ch_vertices_by_pdg)

    unshared_vertices = mg_set ^ ch_set

    # If there's only vertices found that exist in just MadGraph, then
    # create a list of the actual vertices, as instances of MGVertex
    new_vertices = []

    if (unshared_vertices) and not (unshared_vertices & ch_set):
        for i in range(len(mg_verts)):
            v = mg_verts[i]
            parts = v.particles
            if any( part in mg_ghosts for part in parts):
                continue
            vpdg = sorted([ mg_pdg_dict[i] for i in parts ])
            if tuple(vpdg) in unshared_vertices:
                new_vertices.append(v)

    if unshared_vertices:
        if (unshared_vertices & mg_set):
            print("\nThe following vertices, by PDG code, are found in MadGraph but not CalcHEP:")
            print(list(unshared_vertices & mg_set))
        if (unshared_vertices & ch_set):
            print("\nThe following vertices, by PDG code, are found in CalcHEP but not MadGraph:")
            print(list(unshared_vertices & ch_set))
            print("This is not the case ufo2mdl was designed to deal with!")
            print("Please check your input files. However, if you trust the MadGraph files, do:\n")
            print("\tpython ufo_to_mdl <madgraph_folder>\n")
            print("to generate new CalcHEP files from scratch.")
            sys.exit()
        else: 
            print("Writing new CalcHEP files.")
            add_vertices_to_calchep_files(ch_location, new_vertices, param_dict, 
                                          lorentz_dict, part_dict, mg_parts)
    else:
        print("All vertices are consistent...")


"""
WRITING ROUTINES
"""

def write_ch_vertex(vertex, param_dict, lorentz_dict, part_dict, mg_parts, ch_location):
    """
    Writes a CalcHEP vertex, given a MadGraph one.
    """

    # Get the list of particles involved in the vertex
    particles = vertex.particles
    lorentz = vertex.lorentz
    couplings = vertex.couplings

    spin2splus1 = []
    # Get the spins of each particle involved in the vertex
    for part in particles:
        for i in range(len(mg_parts)):
            if mg_parts[i].mgname == part or mg_parts[i].mgantiname == part:
                spin2splus1.append(int(mg_parts[i].spin))
    
    # If there's 4 fermions involved - use the specialist 4-fermion function.
    if spin2splus1.count(2) == 4: 
        return write_four_fermion_vertex(vertex, param_dict, lorentz_dict, 
                                         part_dict, mg_parts, ch_location)
    else:
        return "\n"

    return "\n"

def write_four_fermion_vertex(vertex, param_dict, lorentz_dict, part_dict, mg_parts, ch_location):
    """
    Writes a four-fermion vertex for CalcHEP. This splits the contact interaction into
    two three-body interactions, with a constant propagator. This is achieved by 
    introducing a new auxiliary field.
    """

    """
      - write 4-fermion vertex for calchep using auxiliary field:
 
     ( chibar (gamma) chi ) ( etabar (gamma') eta ) -> entries:
     
     prctls1.mdl
 
     Full name | A     | A+    | PDG   | spinx2 | mass | width | color | aux | tex(A) | tex(A+)
     ...
     <aux>     |~00    |~01    |       |2       |Maux  |0      |1      |!*   | <aux>  | <aux+>
 
     lgrng1.mdl
 
     part 1 | part 2| part 3| part 4| Factor | Lorentz
     ...
     chibar |chi    |~01    |       |i*Maux  |Gamma(m3)
     etabar |eta    |~00    |       |i*Maux  |Gamma'(m3)

    """


    # Need to add new auxiliary particles to the CalcHEP files
    global counter

    newentry =  "{0: <13}".format("aux" + str(counter)) + "|"
    newentry += "{0: <11}".format("~" + str("{:02d}".format(2*counter))) + "|"
    newentry += "{0: <11}".format("~" + str("{:02d}".format(2*counter+1))) + "|"
    newentry +=  "{0: <8}".format("") + "|"
    newentry +=  "{0: <6}".format("2") + "|"
    newentry +=  "{0: <15}".format("Maux") + "|"
    newentry +=  "{0: <15}".format("0") + "|"
    newentry +=  "{0: <5}".format("1") + "|"
    newentry +=  "{0: <3}".format("!*") + "|"
    newentry +=  "{0: <10}".format("aux" + str(counter)) + "|"
    newentry +=  "{0: <14}".format("Aux" + str(counter)) + "|\n"

    print newentry

    new_ch_loc = os.path.abspath(ch_location) + "_ufo2mdl"

    with open(new_ch_loc + "/prtcls1.mdl", 'r') as f, open(new_ch_loc + "/parts_temp", 'w') as g:

        # Copy the original particles
        for line in f:
            g.write(line)

        # Add the new auxiliary particle to the end
        g.write(newentry)

    # Save the new file as the master one
    os.remove(new_ch_loc + "/prtcls1.mdl")
    os.rename(new_ch_loc + "/parts_temp", new_ch_loc + "/prtcls1.mdl")

    # Now write vertices with these auxiliary particles.
    # Create a fake MadGraph particle

    # todo: am here

    parts = mg_parts.copy()
    parts.append(aux_part)

    write_ch_vertex(vertex, param_dict, lorentz_dict, part_dict, parts, ch_location)

    # Increment counter for new particles
    counter+=1

    return "\n"


def write_calchep_files(particles, parameters, vertices, param_dict, lorentz_dict, part_dict, mg_parts):
    """
    Writes CalcHEP files, given all of the information from MadGraph.
    """

    lgrngfile = ""
    varsfile = ""
    prtclsfile = ""
    funcfile = ""
    extlibfile = ""

    for i in range(len(vertices)):
        vertex = vertices[i]
        lgrngfile += write_ch_vertex(vertex, param_dict, lorentz_dict, part_dict, mg_parts, ch_location)

def add_vertices_to_calchep_files(ch_location, new_vertices, param_dict, lorentz_dict, 
                                  part_dict, mg_parts):
    """
    Creates copies of the original CalcHEP files, and adds any new vertices found in the
    MadGraph files to them.
    """

    new_ch_loc = os.path.abspath(ch_location) + "_ufo2mdl"
    if os.path.isdir(new_ch_loc):
        print("\n{} already exists - deleting...").format(new_ch_loc)
        shutil.rmtree(new_ch_loc)        

    shutil.copytree(ch_location, new_ch_loc)

    with open(new_ch_loc+"/lgrng1.mdl", 'r') as f, open(new_ch_loc+"/lgnrg_temp", 'w') as g:

        # Copy over the vertices we want to keep
        for line in f:
            g.write(line)

        # Now write all the new vertices in.
        for i in range(len(new_vertices)):
            vertex = new_vertices[i]
            g.write(write_ch_vertex(vertex, param_dict, lorentz_dict, part_dict, mg_parts, ch_location))


    print("New CalcHEP vertices added to {}.").format(new_ch_loc)

"""
Main functions
"""

def convert(mg_location):
    """
    Create CalcHEP files from a given set of MadGraph files.
    """
    
    # Check the MadGraph files are where the user says
    check_mg_files(mg_location)

    print("Reading MadGraph files...")

    # Read in the MadGraph files
    params = import_mg_params(mg_location)
    parts = import_mg_parts(mg_location)
    couplings = import_mg_couplings(mg_location)
    lorentz = import_mg_lorentz(mg_location)
    verts = import_mg_verts(mg_location)

    print("Done.")

    print("Writing CalcHEP output...")

    write_calchep_files(ch_location, particles, parameters, vertices)

    print("Done.")
    
    
def compare(mg_location, ch_location):
    """
    Compare a list of vertices from MadGraph with those from CalcHEP.
    By this, we simply check that each vertex is present in each set of 
    model files.
    """
    
    # Check both sets of files are where they should be
    check_mg_files(mg_location)
    check_ch_files(ch_location)

    print("Reading MadGraph files...")
    
    # Read in the MadGraph files
    mg_params = import_mg_params(mg_location)
    mg_parts = import_mg_parts(mg_location)
    mg_couplings = import_mg_couplings(mg_location)
    mg_lorentz = import_mg_lorentz(mg_location)
    mg_verts = import_mg_verts(mg_location)

    print("Reading CalcHEP files...")

    # Read in the CalcHEP files
    ch_params = import_ch_params(ch_location)
    ch_parts = import_ch_parts(ch_location)
    ch_verts = import_ch_verts(ch_location)

    print("Done.")
    print("Comparing MadGraph and CalcHEP files...")
    print("Checking external variables...")

    # CalcHEP can compute widths 'on the fly' - these will be parameters in MadGraph, so add
    # these as another group...
    ch_widths = [x.width.strip('!') for x in ch_parts if x.width.startswith('!')]

    shared_by_name = set([x.name for x in ch_params] + ch_widths) & set([x.name for x in mg_params])
    unshared_by_name = set([x.name for x in ch_params] + ch_widths) ^ set([x.name for x in mg_params])

    # Go through those with unmatching names and add the *full* parameter member instance to it
    unshared_full = []
    for i in range(len(mg_params)):
        if mg_params[i].name in unshared_by_name:
            unshared_full.append(mg_params[i])
    for i in range(len(ch_params)):
        if ch_params[i].name in unshared_by_name:
            unshared_full.append(ch_params[i])

    # For those unshared parameters, run them through the dictionary of common misnomers to 
    # check they're actually not repeated
    param_dict = get_mg_ch_param_dict()
    mg_full, ch_full = run_through_parameter_dict(unshared_full, param_dict)

    # If anything is found in MadGraph, tell the user
    if len(mg_full) > 0:
        print("The following external parameters are defined in MadGraph but not CalcHEP:")
        print([x.name() for x in mg_full])
        if len(ch_full) == 0: sys.exit()
    # If anything is found in CalcHEP, tell the user, and kill.
    if len(ch_full) > 0:
        print("The following external parameters are defined in CalcHEP but not MadGraph:")
        print([x.name() for x in mg_full])
        sys.exit()

    print("All external variables are consistent!")
    print("Checking particle content...")

    # Obtain a dictionary to translate between MG and CH particle names, just in case they are
    # different. They probably won't be.
    part_dict = get_mg_ch_parts_dict(mg_parts, ch_parts)

    print("All particles are consistent...")
    print("Checking vertices...")

    lorentz_dict = get_mg_ch_lorentz_dict()

    compare_vertices(mg_verts, mg_parts, mg_params,
                     ch_verts, ch_parts, ch_params,
                     ch_location, param_dict, lorentz_dict, part_dict)

    print("Done!")
    
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

    print("****************************************************************")
    print("Welcome to the .ufo to .mdl file converter!")
    print("If you use this, please cite the GUM manual: arXiv:19XX.YYYYY.")
    print("****************************************************************\n")
    
    if len(sys.argv) == 1 or len(sys.argv) > 3:
        usage()
        
    elif len(sys.argv) == 2:
        convert(sys.argv[1])
        
    elif len(sys.argv) == 3:
        compare(sys.argv[1], sys.argv[2])

"""
 TODO
 
  - add a --force flag so if the files don't agree, automatically produce CH files from MG
  - write 4-fermion vertex for calchep using auxiliary field:
 
     ( chibar (gamma) chi ) ( etabar (gamma') eta ) -> entries:
     
     prctls1.mdl
 
     Full name | A     | A+    | PDG   | spinx2 | mass | width | color | aux | tex(A) | tex(A+)
     ...
     <aux>     |~00    |~01    |       |2       |Maux  |0      |1      |!*   | <aux>  | <aux+>
 
     lgrng1.mdl
 
     part 1 | part 2| part 3| part 4| Factor | Lorentz
     ...
     chibar |chi    |~01    |       |i*Maux  |Gamma(m3)
     etabar |eta    |~00    |       |i*Maux  |Gamma'(m3)
 
 
 """