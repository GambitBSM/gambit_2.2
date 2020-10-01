#  GUM: GAMBIT Universal Models
#  ****************************
#  \file
#
#  Routines for writing CalcHEP (.mdl) files from MadGraph (.ufo) files.
#
#  *************************************
#
#  \author Sanjay Bloor
#          (sanjay.bloor12@imperial.ac.uk)
#  \date 2018, 2019
#
#  **************************************

"""
The user can either compare between .ufo and .mdl files to check everything
is consistent between all files:
    
    python ufo_to_mdl <madgraph_folder> <calchep_folder>
    
Or just run on a set of .ufo files and generate .mdl files from scratch:

    python ufo_to_mdl <madgraph_folder>
   
If you use this script as a standalone please cite the GUM manual:
      arXiv:20XX.YYYYY

"""

import sys
import os
import re
import shutil
from collections import Counter
from future.utils import iteritems

counter = 0
opt = 0
optimisations = {}

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
        self.spin = spin # 2s+1
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
    Internal object representing a Lorentz structure
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
    coups = re.findall(r'(.*?) = Coupling\((.*?)\)\n', lines)

    for i in range(len(coups)):

        mgname = coups[i][0]
        
        # Match 'value = VALUE'
        value = re.findall(r'value = \'(.*?)\'', str(coups[i]))[0]
        order = re.findall(r'order = {(.*?)}', str(coups[i]))[0]

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
        # Remove whitespace, quote marks and the preceding "L."
        lorentz = [x.strip().strip('L.') for x in lorentz]

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


def translate_lorentz_structure(lorentz):
    """
    Translates the Lorentz structure of a MadGraph vertex into the language 
    of CalcHEP.
    """

    # todo: factor of i for goldstones? can't have complex numbers in CH.
    # todo: epsilon(a,b,c,d)
    # todo: minus signs in gammas..?


    # If it's a 4-fermion vertex, we're going to have to split it up.
    if lorentz.mgname[:4] == "FFFF":
        return translate_four_fermion_lorentz(lorentz)

    ch_struct = lorentz.structure

    # Projection operators in MadGraph are ProjP(a,b) = (I + gamma^5)_{a,b}/2, (M<->P, +<->-)
    ch_struct = re.sub(r'ProjP\((.*?),(.*?)\)', 'half*(1+G5)', ch_struct)
    ch_struct = re.sub(r'ProjM\((.*?),(.*?)\)', 'half*(1-G5)', ch_struct)

    # Momenta in MadGraph are P(A,N) = P_N^A
    # Momenta in CalcHEP are pN.mA
    ch_struct = re.sub(r'P\((.*?),(.*?)\)', r'p\2.m\1', ch_struct)

    # Metric in MadGraph is Metric(mu,nu) = g_{mu,nu}
    # Metric in CalcHEP is m1.m2 = g_{m1,m2}
    ch_struct = re.sub(r'Metric\((.*?),(.*?)\)',r'm\1.m\2', ch_struct)

    # Dirac matrix in MadGraph is Gamma(mu,a,b) = gamma^mu_{a,b}
    # Dirac matrix in CalcHEP is G(m1) -> gamma^m1
    ch_struct = re.sub(r'Gamma\((.*?),(.*?),(.*?)\)',r'G(m\1)', ch_struct)

    # Gamma5 in MadGraph is Gamma5(a,b) = gamma^5_{a,b}
    # Gamma5 in CalcHEP is G5
    ch_struct = re.sub(r'Gamma5\((.*?),(.*?)\)','G5', ch_struct)

    # Now we have to split this all up into Dirac algebra and the other stuff.

    # TODO: Do some algabra on gamma matrices? Clean it up? Leave it as one big string?
    return ch_struct

def translate_four_fermion_lorentz(lorentz):
    """
    Translates the Lorentz structure of a MadGraph vertex into the language 
    of CalcHEP, for a 4-fermion vertex.
    """

    ch_struct = lorentz.structure

    # These vertices are like:
    # ~ (chibar gamma chi) (etabar gamma' eta)
    # To be split into:
    # chibar (gamma) chi AUX   +   etabar gamma' eta (AUXbar)
    # to be denoted left + right
    # The vertices will have already been organised s.t. only those 
    # with consistent Lorentz indices are combined

    # Use the multiplication to split the string up
    strings = ch_struct.split('*')

    # These are now strings like Gamma5(2,1) or Gamma(-1,4,3)
    # Find all of those with both a +1 and +2 in, and same for +3 and +4
    # Those with a +1 and a +2 are on the left hand side, for sure.
    # Those with a +3 and a +4 are on the right hand side.
    left = []
    right = []

    # This procedure works for 4-fermion Lorentz structures with 
    # 4 or fewer Lorentz indices (i.e. if there are no tensor currents)
    if len(strings) <= 4:
        for i in strings:
            # Grab the indices labelling the gamma matrices
            nums = re.search(r'\((.*?)\)',i).group(1).split(',')
            # Make these integers
            nums = [int(x) for x in nums]
            if 1 in nums or 2 in nums:
                left.append(i)
            elif 3 in nums or 4 in nums:
                right.append(i)
    # If there's 5 different structures....
    elif len(strings) == 5:
        # We probably have contractions with a gamma5 in this case.
        # Just check explicitly:
        search = [r"Gamma\(", r"Gamma5\("]
        to_search = re.compile('|'.join(sorted(search, key=len, reverse=True)))
        matches = (to_search.search(el) for el in strings)
        counts = Counter(match.group() for match in matches if match)
                
        # Bingo! 
        if counts["Gamma("] == 4 and counts["Gamma5("] == 1:
            # Now need to figure out which order the Lorentz indices go.
            # Doesn't matter which side has the g5 - it commutes with Sigma

            # Find the gamma matrix attached to the fermions labelled 1 and 4
            first_index = 0
            last_index = 0

            for i in strings:
                # Grab the indices labelling the gamma matrices
                nums = re.search(r'\((.*?)\)',i).group(1).split(',')
                # Make these integers
                nums = [int(x) for x in nums]

                if 2 in nums:
                    first_index = nums[0]
                if 3 in nums:
                    last_index = nums[0]

            # Just put the CalcHEP syntax in here.
            if first_index == 0 or last_index == 0:
                print("Could not find the correct indices.")
                print("Please check your UFO files.")
                sys.exit()
            # (mu,nu)(nu,mu)
            if first_index == last_index:
                left.append("G5*G(m3)*G(M3)")
                right.append("G(M3)*G(m3)")
            # (mu,nu)(mu,nu)
            else:
                left.append("G5*G(m3)*G(M3)")
                right.append("G(m3)*G(M3)")

        # *Pretty* sure this is the only dim 5 operator we'll get...
        else:
            print("Unexpected Lorentz structure found.")
            print("Quitting...")
            sys.exit()

    # Or it's just something super weird I've not come across either.
    # Don't think we'll ever get dimension 6, can't think of how.
    # Note that if you have (chi g5 sigma chi)(eta g5 sigma eta) this
    # simply equals (chi sigma chi) (eta sigma eta): go rewrite your vertices!
    else:
        print("Unexpected Lorentz structure found.")
        print("Quitting...")
        sys.exit()


    # Now want to perform operations on left and right
    left = '*'.join(left)
    right = '*'.join(right)
    
    # Now we do magic on the Lorentz structures.

    # Dirac matrix in MadGraph is Gamma(mu,a,b) = gamma^mu_{a,b}
    # Dirac matrix in CalcHEP is G(m1) -> gamma^m1

    # This always gets put on the auxiliary particle so it carries the
    # correct Lorentz structure in the propagator...
    
    # Let's do the case where we have two Gammas first, otherwise
    # the single routines will give incorrect indices

    # Order of Lorentz indices is important, so make sure we follow them through
    if ( re.search(r'Gamma\((.*?),(.*?),(.*?)\)*Gamma\((.*?),(.*?),(.*?)\)', left) and
         re.search(r'Gamma\((.*?),(.*?),(.*?)\)*Gamma\((.*?),(.*?),(.*?)\)', right) ):
            
        first_index = 0
        last_index = 0

        nums_l = re.findall(r'\((.*?)\)',left)
        nums_r = re.findall(r'\((.*?)\)',right)

        n_l = []
        n_r = []
        for i in nums_l:
            n_l.append([int(x) for x in i.split(',')])
        for i in nums_r:
            n_r.append([int(x) for x in i.split(',')])

        for i in n_l:
            if 2 in i:
                first_index = i[0]
        for j in n_r:
            if 3 in j:
                last_index = j[0]

        if first_index == 0 or last_index == 0:
            print("Could not find the correct indices.")
            print("Please check your UFO files.")
            sys.exit()
        # (mu,nu)(nu,mu)
        if first_index == last_index:
            left = re.sub(r'Gamma\((.*?),(.*?),(.*?)\)*Gamma\((.*?),(.*?),(.*?)\)',r'G(m3)*G(M3)', left)
            right = re.sub(r'Gamma\((.*?),(.*?),(.*?)\)*Gamma\((.*?),(.*?),(.*?)\)',r'G(M3)*G(m3)', right)
        # (mu,nu)(mu,nu)
        else:
            left = re.sub(r'Gamma\((.*?),(.*?),(.*?)\)*Gamma\((.*?),(.*?),(.*?)\)',r'G(m3)*G(M3)', left)
            right = re.sub(r'Gamma\((.*?),(.*?),(.*?)\)*Gamma\((.*?),(.*?),(.*?)\)',r'G(m3)*G(M3)', right)

    
    # left = re.sub(r'Gamma\((.*?),(.*?),(.*?)\)*Gamma\((.*?),(.*?),(.*?)\)',r'G(m3)*G(M3)', left)
    # right = re.sub(r'Gamma\((.*?),(.*?),(.*?)\)*Gamma\((.*?),(.*?),(.*?)\)',r'G(m3)*G(M3)', right)

    # Now the case where there's one. 
    left = re.sub(r'Gamma\((.*?),(.*?),(.*?)\)',r'G(m3)', left)
    right = re.sub(r'Gamma\((.*?),(.*?),(.*?)\)',r'G(m3)', right)

    # Gamma5 in MadGraph is Gamma5(a,b) = gamma^5_{a,b}
    # Gamma5 in CalcHEP is G5
    left = re.sub(r'Gamma5\((.*?),(.*?)\)','G5', left)
    right = re.sub(r'Gamma5\((.*?),(.*?)\)','G5', right)

    # Identity is 1!
    left = re.sub(r'Identity\((.*?),(.*?)\)','1', left)
    right = re.sub(r'Identity\((.*?),(.*?)\)','1', right)

    # Sigma matrices don't seem to be defined in CalcHEP, so write them out by hand
    # todo: factor of i needs to be added to the vertex factor!
    left = re.sub(r'Sigma\((.*?),(.*?)\)','(half*(G(m3)*G(M3)-G(M3)*G(m3)))', left)
    right = re.sub(r'Sigma\((.*?),(.*?)\)','(half*(G(m3)*G(M3)-G(M3)*G(m3)))', right)

    return left, right


def get_mg_ch_lorentz_dict(mg_lorentz):
    """
    Returns a dictionary to translate between Lorentz structures in MadGraph
    and CalcHEP.
    """

    ld = {}

    for i in range(len(mg_lorentz)):
        ld[mg_lorentz[i].mgname] = translate_lorentz_structure(mg_lorentz[i])

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
        if not mg_parts[i].selfconj: mg_pairs[-1*int(mg_parts[i].pdg_code)] = mg_parts[i].antiname

    keys = []
    values = []

    for i in range(len(ch_parts)):
        keys.append( mg_pairs.get(ch_parts[i].pdg_code) )
        values.append( ch_parts[i].name )
        if not ch_parts[i].selfconj: 
            keys.append( mg_pairs.get(-1*int(ch_parts[i].pdg_code) ))
            values.append( ch_parts[i].antiname)

    # Now construct a dictionary of {'MadGraph name': 'CalcHEP name'}
    part_dict = dict(zip(keys, values))

    return part_dict

def get_from_part_dict(particle_name, part_dict):
    """
    Returns the CalcHEP particle from the MadGraph particle.

    If we are adding a new 4-fermion interaction there will be an auxiliary
    particle which does not exist in the MadGraph files, so won't be present
    in the dictionary. So let's deal with this.
    """

    if re.match(r'~\d+', particle_name):
        return particle_name
    else:
        return part_dict.get(particle_name)

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

def compare_vertices(mg_verts, mg_parts, mg_params,
                     mg_lorentz, mg_couplings,
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

        # Find ghosts and goldstones and tensor colour structures
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

    # 4-gluon vertex -- hard-coded as it is in SARAH/FeynRules.
    # This interaction is split into 2 3-body interactions with a coloured 
    # tensor auxiliary propagator. Let's see if that is in the CH files.
    if (21, 21, 21, 21) in unshared_vertices:

        # Let's get the gluon name
        gluon = ""
        for i in range(len(ch_parts)):
            p = ch_parts[i]
            if int(p.pdg_code) == 21:
                gluon = p.name
        
        # We want [g, g, g.t] vertex...
        lookup = tuple([gluon, gluon, gluon+'.t'])
        for i in range(len(ch_verts)):
            v = ch_verts[i]
            parts = v.particles
            # If we find what we want, remove the 4G vertex
            if lookup == tuple(sorted(parts)):
                unshared_vertices.remove( (21, 21, 21, 21) )
                break

        # If it's not found, then we can add it later.


    # If there's only vertices found that exist in just MadGraph, then
    # create a list of the actual vertex objects (i.e. instances of MGVertex)
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

    # 
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
                                          lorentz_dict, part_dict, mg_parts, 
                                          mg_lorentz, mg_couplings)
            return False
    else:
        print("All vertices are consistent...")
        return True


"""
WRITING ROUTINES
"""

def get_lorentz_index(lorentz_structure, vertex_couplings, index, mg_lorentz):
    """
    Finds out the Lorentz index of a given interaction. 
    """

    for i in range(len(mg_lorentz)):
        l = mg_lorentz[i]

        # We've found the right structure according to the internal name.
        if l.mgname == lorentz_structure:

            # Find the corresponding entry in the vertex coupling dict.
            tomatch = '(0,{})'.format(index)

            coup = ""
            for k in vertex_couplings:
                if tomatch in k:
                    coup = k
            if not coup:
                print("Coupling not found. Aborting.")
                sys.exit()

            # Count the number of Lorentz indices.
            numgammas = l.structure.split(' ')[0].count("Gamma(")
            numsigmas = l.structure.split(' ')[0].count("Sigma(")

            # No Lorentz indices
            if numgammas + numsigmas == 0:
                return 0

            elif numgammas == 2 and numsigmas == 0:
                return 1

            elif numgammas == 4 and numsigmas == 0:
                return 2

            elif numgammas == 0 and numsigmas == 2:
                return 2

            else: 
                print("Weird number of Lorentz indices... I'm out!")
                sys.exit()



def write_ch_vertex(vertex, param_dict, lorentz_dict, part_dict, mg_parts, 
                    ch_location, mg_lorentz, mg_couplings):
    """
    Writes a CalcHEP vertex, given a MadGraph one.
    """

    new_ch_loc = os.path.abspath(ch_location) + "_ufo2mdl"

    # Get the list of particles involved in the vertex
    particles = vertex.particles
    lorentz = vertex.lorentz
    couplings = vertex.couplings

    spin2splus1 = []
    pdgs = []
    parts = []

    # Get the spins of each particle involved in the vertex
    for part in particles:
        SET = False # Has the particle been added?
        for i in range(len(mg_parts)):
            if mg_parts[i].mgname == part:
                spin2splus1.append(int(mg_parts[i].spin))
                pdgs.append(int(mg_parts[i].pdg_code))
                parts.append(mg_parts[i].name)
                SET = True
            elif mg_parts[i].mgantiname == part:
                pdgs.append(-1*int(mg_parts[i].pdg_code))
                spin2splus1.append(int(mg_parts[i].spin))
                parts.append(mg_parts[i].antiname)
                SET = True
        ## 4-fermion routines will give an auxiliary particle called
        ## something like ~0, ~1, etc. Add these even though they're not
        ## in the MadGraph files, otherwise this will break
        if SET == False and re.match(r'~\d+', part):
            parts.append(part)

    # If there's 4 fermions involved - use the specialist 4-fermion function.
    if spin2splus1.count(2) == 4:
        write_four_fermion_vertex(vertex, param_dict, lorentz_dict, 
                                  part_dict, mg_parts, mg_lorentz, 
                                  mg_couplings, ch_location)

            
    # If it's a 4-gluon vertex, this is a hard-coded result: 
    # just use the specialist function
    elif pdgs == [21, 21, 21, 21]:
        return write_four_gluon_vertex(part_dict, mg_parts, ch_location)
    
    #todo: For nice formatting, get the numbers of spaces for each header.
    # Otherwise - let's go!
    else:
        coup = ""
        lor = ""
        # If we are adding a 4-fermion vertex, the conversion to CalcHEP-friendly language will
        # have already been done, so no need to do it again.
        conversion_done = False
        for part in particles:
            if re.match(r'~\d+', part):
                conversion_done = True

        if conversion_done:
            lor = lorentz
            coup = couplings
        else:
            # Rewrite the Lorentz structure, and save each coupling
            lor = []
            coups = []
            for i in range(len(mg_lorentz)):
                l = mg_lorentz[i]
                for j in range(len(vertex.lorentz)):
                    if vertex.lorentz[j] == l.mgname:
                        s = mg_lorentz[i].mgname
                        # Rewrite the Lorentz structure into CalcHEP language
                        lor.append(lorentz_dict.get(s))
                        coups.append(vertex.couplings[j])

            # Rewrite the couplings in CH-friendly format

            # Go through the list of couplings, and swap in their actual values
            coups_dict = {}
            for i in coups:

                # The key will look like (i,j) -> i-th entry of the colour tensor space, 
                # j-th entry of the lorentz tensor space. Since we're only dealing with 
                # cases where we have the identity in colour space (for now), we can ditch
                # the first entry of the key, and just keep the jth value.
                key, value = i.split(':')[0].split(',')[1][:-1] , i.split(':')[1][2:]
                coups_dict[ key ] = value

            # Now we need to pair up each element of the lorentz array
            # with the correct coupling.

            lc_dict = {}
            for k, v in iteritems(coups_dict):
                for j in mg_couplings:
                    if v == j.mgname:
                        lc_dict[ lor[int(k)] ] = j.value

            # Dictionary of lorentz structures & couplings -- updated to be in C 
            # syntax, and only the ones relevant for the interactions at hand.
            c_dict = c_ify_couplings(mg_couplings, param_dict)

            lc_dict_updated = {}
            
            for k, v in iteritems(lc_dict):
                lc_dict_updated.update({k: c_dict[v]})

            # Put all the couplings in the Lorentz bit as one big horrible string
            horrible_string = []
            for k,v in iteritems(lc_dict_updated):
                horrible_string.append( v+"*"+k[0] ) # Coupling * Lorentz factor

            # CalcHEP doesn't understand +-, so we
            # can't simply use the string.join(list) method,
            # as the first character might be a "-"
            lor = horrible_string[0]
            for i in horrible_string[1:]:
                if i.startswith("-"): 
                    lor += i
                else:
                    lor += "+"+i
            coup = "i" # TODO should this be i or 1?
        
        # Now we can do some writing.
        vertex = ""
        vertex += "{0: <13}".format(get_from_part_dict(parts[0], part_dict)) + "|"
        vertex += "{0: <13}".format(get_from_part_dict(parts[1], part_dict)) + "|"
        vertex += "{0: <13}".format(get_from_part_dict(parts[2], part_dict)) + "|"
        if len(parts) == 4:
            vertex += "{0: <13}".format(get_from_part_dict(parts[3], part_dict)) + "|"
        else:
            vertex += "{0: <13}".format("") + "|"

        vertex += "{0: <31}".format(coup) + "|"
        vertex += "{}\n".format(lor)

        with open(new_ch_loc + "/lgrng1.mdl", 'r') as f, open(new_ch_loc + "/lgrng_temp", 'w') as g:

            # Copy the original vertices
            for line in f:
                g.write(line)

            # Add the new vertex to the end
            g.write(vertex)

        # Save the new file as the master one
        os.remove(new_ch_loc + "/lgrng1.mdl")
        os.rename(new_ch_loc + "/lgrng_temp", new_ch_loc + "/lgrng1.mdl")

        return vertex

def c_ify_couplings(mg_couplings, param_dict):
    """
    Returns a dictionary mapping Python syntax to C syntax for the couplings.
    """

    coupling_dict = {}
    # Go through the couplings and make the MG Python syntax into CH C syntax
    for i in mg_couplings:
        

        # create a local copy to play about with
        t = i.value

        # cmath.sqrt(2) --> Sqrt2 
        t = re.sub(r'cmath.sqrt\(2\)', 'Sqrt2', t)
        # complex(0,1) --> i
        t = re.sub(r'complex\(0,1\)', 'i', t)
        # cmath.pi --> Pi
        t = re.sub(r'cmath.pi', 'Pi', t)
        # 1/a**b --> 1/a^b
        t = re.sub(r'/(\w+)\*\*(\w+)', r'/\1^\2', t)
        # a**b --> a^b
        t = re.sub(r'(\w+)\*\*(\w+)', r'\1^\2', t)
        # Any floating point number with a dot after it -> remove the dot
        # e.g. 2.*Lambda -> 2*Lambda. Obviously not if it's a number like 1.5.
        t = re.sub(r'(\d+)\.([^0-9]+)', r'\1\2', t)

        # ComplexConjugate
        # TODO - since CalcHEP can't hack complex parameters, will want to 
        # either:
        # a) split the complex parameter into real and imaginary parts;
        # b) throw an error;
        # c) assume everything is okay (bad idea -- but doing this for now);
        t = re.sub(r'complexconjugate\((\w+)\)', r'\1', t)

        # Finally go through the param_dict and change any parameters
        # that are named differently between MadGraph and CalcHEP
        for k, v in iteritems(param_dict):
            t = re.sub(k, v, t)
        
        # Look, another dictionary!
        coupling_dict[i.value] = t

    return coupling_dict

def write_four_fermion_vertex(vertex, param_dict, lorentz_dict, part_dict,
                              mg_parts, mg_lorentz, mg_couplings, ch_location):
    """
    Writes a four-fermion vertex for CalcHEP. This splits the contact 
    interaction into two three-body interactions, with a constant propagator. 
    This is achieved by introducing a new auxiliary field.

    A vertex from MadGraph can have multiple entries, for instance the 
    Z-f-f vertex from the Standard Model would have a separate entry for
    the vector and axial vector interactions. We'll split these up in CalcHEP,
    so each contact interaction corresponds to one of these Dirac operators.
    
    ///////////////////////////////////////////////////////////////////////////
    
    There are 3 categories to split these interactions up into: 
    

    %%% CASE 1 %%%
    
    no lorentz indices --> scalar mediator

    e.g. 
       ( chibar chi ) ( etabar eta )

    would have the following entries:

    prctls1.mdl
    Full name | A     | A+    | PDG   | spinx2 | mass | width | color | aux | tex(A) | tex(A+)
    ..
    <aux>     |~00    |~01    |       |0       |Maux  |0      |1      |!*   | <aux>  | <aux+>
    
    lgrng1.mdl
    part 1 | part 2| part 3| part 4| Factor | Lorentz
    ...
    chibar |chi    |~01    |       |Maux    |1
    etabar |eta    |~00    |       |Maux    |1

    %%% CASE 2 %%% 
    
    1 lorentz index --> vector mediator

    e.g. 
        chibar (gamma^mu) chi ) ( etabar (gamma'_mu) eta )

    would have the following entries:

    prctls1.mdl
    Full name | A     | A+    | PDG   | spinx2 | mass | width | color | aux | tex(A) | tex(A+)
    ..
    <aux>     |~00    |~01    |       |2       |Maux  |0      |1      |!*   | <aux>  | <aux+>
    
    lgrng1.mdl
    part 1 | part 2| part 3| part 4| Factor | Lorentz
    ...
    chibar |chi    |~01    |       |i*Maux  |Gamma(m3)
    etabar |eta    |~00    |       |i*Maux  |Gamma'(m3)

    %%% CASE 3. %%%

    2 lorentz indices --> tensor mediator
    
    e.g. 
        (chibar gamma^mu gamma^nu chi) (etabar gamma_mu gamma_nu eta)

    would have the following entries:

    prctls1.mdl
    Full name | A     | A+    | PDG   | spinx2 | mass | width | color | aux | tex(A) | tex(A+)
    ..
    <aux>     |~00    |~01    |       |2       |Maux  |0      |1      |!*   | <aux>  | <aux+>
    
    lgrng1.mdl
    part 1 | part 2| part 3| part 4| Factor | Lorentz
    ...
    chibar |chi    |~01.t  |       |1         |G(m3)*G(M3)
    etabar |eta    |~00.t  |       |1         |G(m3)*G(M3)
    """

    # Need to add new auxiliary particles to the CalcHEP files
    global counter

    # For some optimisations to save
    global opt
    global optimisations

    # Location of new files
    new_ch_loc = os.path.abspath(ch_location) + "_ufo2mdl"

    # TODO: more complex colour structures
    # For now -- if the colour structure is non-trivial, UFO2MDL will complain
    for i in vertex.color[:]:
        if i not in ["1", "Identity(3,4)"]:
            print("ufo2mdl does not yet support funky color structures. Sorry.")
            sys.exit()

    # Number of new vertices we are adding -- the number of operators we
    # have for each vertex; a 'subvertex' if you like.
    numvertices = len(vertex.lorentz)

    # CalcHEP-ify the particle names - same for each sub-vertex of course
    v = []
    for i in vertex.particles[:]:
        for j in range(len(mg_parts)):
            if mg_parts[j].mgname == i:
                v.append( part_dict.get(mg_parts[j].name) )
            elif mg_parts[j].mgantiname == i:
                v.append( part_dict.get(mg_parts[j].antiname) )

    # Construct a dictionary of all coupling - lorentz pairs
    # This looks like ( key : value ) = ( MG lorentz struct : MG coupling )
    lorentz_couplings = {}

    for i in range(numvertices):
        # find the coupling with the index (0, i)
        tomatch = "(0,{})".format(i)
        for coupling in vertex.couplings:
            if tomatch in coupling:
                v = vertex.lorentz[i]
                entry = coupling.split(':')[1][2:]
                lorentz_couplings[v] = entry

    # Now we have a dictionary of each pair of (left, right) fermion operators
    # and the associated coupling. 

    # Now get a dictionary of lorentz structures & couplings -- updated to be 
    # in C syntax, and only the ones relevant for the interactions at hand.
    c_dict = c_ify_couplings(mg_couplings, param_dict)
    
    # Now go through the MG couplings and create a dictionary of 
    # (key : value ) = ( MG coupling index : C expression for the coupling ) 
    coupling_dict = {}
    for i in mg_couplings:
        coupling_dict[i.mgname] = c_dict[i.value]

    # Save each coupling as a global optimised function
    for k2, v2 in iteritems(lorentz_couplings):
        newstring = "GUM{:0>3}".format(str(opt))
        opt = opt+1
        optimisations[ newstring ] = coupling_dict[v2]

    # for k1, v1 in iteritems(coupling_dict):
    #     #print k1, v1

    #     newstring = "GUM{:0>3}".format(str(opt))
    #     opt = opt+1
    #     optimisations[ newstring ] = v1

    # Got all we need!
    # Go through each sub-vertex and add it individually
    for i in range(numvertices):

        # Find out how many Lorentz indices we're working with, so the 
        # mediator particle carries the right current
        lorentz_indices = get_lorentz_index(vertex.lorentz[i], vertex.couplings,
                                            i, mg_lorentz)

        # Lorentz structure as known by MG
        vli = vertex.lorentz[i]
        # Tuple of Lorentz structure for each auxiliary interaction
        lorentz = lorentz_dict[vli]
        # Corresponding coupling in C form
        coupling = coupling_dict[lorentz_couplings[vli]]

        # Get the optimised name for the coupling
        optimisedout = ""
        for k, v in iteritems(optimisations):
            if v == coupling: 
                optimisedout = v + "*Maux"
        # todo error

        newentry = ""

        # Add a different particle based on the number of Lorentz indices
        if lorentz_indices == 0:
            
            newentry += "{0: <13}".format("aux" + str(counter)) + "|"
            newentry += "{0: <11}".format("~" + str("{:02d}".format(2*counter))) + "|"
            newentry += "{0: <11}".format("~" + str("{:02d}".format(2*counter+1))) + "|"
            newentry += "{0: <8}".format("") + "|"
            newentry += "{0: <6}".format("0") + "|"
            newentry += "{0: <15}".format("Maux") + "|"
            newentry += "{0: <15}".format("0") + "|"
            newentry += "{0: <5}".format("1") + "|"
            newentry += "{0: <3}".format("!*") + "|"
            newentry += "{0: <10}".format("aux" + str(counter)) + "|"
            newentry += "{0: <14}".format("Aux" + str(counter)) + "|\n"

        elif lorentz_indices == 1:

            newentry += "{0: <13}".format("aux" + str(counter)) + "|"
            newentry += "{0: <11}".format("~" + str("{:02d}".format(2*counter))) + "|"
            newentry += "{0: <11}".format("~" + str("{:02d}".format(2*counter+1))) + "|"
            newentry += "{0: <8}".format("") + "|"
            newentry += "{0: <6}".format("2") + "|"
            newentry += "{0: <15}".format("Maux") + "|"
            newentry += "{0: <15}".format("0") + "|"
            newentry += "{0: <5}".format("1") + "|"
            newentry += "{0: <3}".format("!*") + "|"
            newentry += "{0: <10}".format("aux" + str(counter)) + "|"
            newentry += "{0: <14}".format("Aux" + str(counter)) + "|\n"

        elif lorentz_indices == 2:

            # First particle (real?)
            newentry += "{0: <13}".format("aux" + str(counter)) + "|"
            newentry += "{0: <11}".format("~" + str("{:02d}".format(2*counter))) + "|"
            newentry += "{0: <11}".format("~" + str("{:02d}".format(2*counter+1))) + "|"
            newentry += "{0: <8}".format("") + "|"
            newentry += "{0: <6}".format("2") + "|"
            newentry += "{0: <15}".format("Maux") + "|"
            newentry += "{0: <15}".format("0") + "|"
            newentry += "{0: <5}".format("1") + "|"
            newentry += "{0: <3}".format("!*") + "|"
            newentry += "{0: <10}".format("aux" + str(counter)) + "|"
            newentry += "{0: <14}".format("Aux" + str(counter)) + "|\n"

        else:

            print(("UFO2MDL does not know how to work with {} lorentz "
                   "indices.\n").format(lorentz_indices))
            sys.exit()


        with open(new_ch_loc + "/prtcls1.mdl", 'r') as f, open(new_ch_loc + "/parts_temp", 'w') as g:

            # Copy the original particles
            for line in f:
                g.write(line)

            # Add the new auxiliary particle to the end
            g.write(newentry)

        # Save the new file as the master one
        os.remove(new_ch_loc + "/prtcls1.mdl")
        os.rename(new_ch_loc + "/parts_temp", new_ch_loc + "/prtcls1.mdl")

        # By convention, assign the coupling to the first vertex, then just give
        # the other vertex the correct factor given the Lorentz structure.

        # Let's create a new MG vertex object and use the existing function
        # to create a CH vertex
        if lorentz_indices in [0, 1]:
            v1parts = [vertex.particles[0], vertex.particles[1], "~" + str("{:02d}".format(2*counter))]
            v2parts = [vertex.particles[2], vertex.particles[3], "~" + str("{:02d}".format(2*counter+1))]
        # Need to indicate that the auxiliary particle has tensor indices
        elif lorentz_indices == 2:
            v1parts = [vertex.particles[0], vertex.particles[1], "~" + str("{:02d}.t".format(2*counter))]
            v2parts = [vertex.particles[2], vertex.particles[3], "~" + str("{:02d}.t".format(2*counter+1))]

        # TODO check this!!!
        if lorentz_indices == 0:
            v2coups = "-i*Maux"
        elif lorentz_indices == 1:
            v2coups = "-i*Maux"
        elif lorentz_indices == 2:
            v2coups = "i"

        vertex1 = MGVertex("vertex1", v1parts, "", lorentz[0], optimisedout)
        vertex2 = MGVertex("vertex2", v2parts, "", lorentz[1], v2coups)

        # Try again with the two 3-body vertices
        write_ch_vertex(vertex1, param_dict, lorentz_dict, part_dict, mg_parts, 
                        ch_location, mg_lorentz, mg_couplings)
        write_ch_vertex(vertex2, param_dict, lorentz_dict, part_dict, mg_parts, 
                        ch_location, mg_lorentz, mg_couplings)

        # Increment counter for auxiliary particles
        counter+=1

    return "\n"

def write_four_gluon_vertex(part_dict, mg_parts, ch_location):
    """
    Writes the Standard Model 4-gluon interaction in CalcHEP via an auxiliary
    tensor propagator carrying colour charge.  

    # todo
    """

    # Get the name of the gluon in CalcHEP

    # v =[]
    # for i in mg_parts: if mgp[i].pdg == 21: glu=mgp[i].mgname
    # chname = part_dict.get(glu)
    # vertex = [chname, chname, chname+'.t']
    # etc.

    return "\n"

def write_calchep_files(particles, parameters, vertices, param_dict, 
                        lorentz_dict, part_dict):
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
        lgrngfile += write_ch_vertex(vertex, param_dict, lorentz_dict, part_dict, particles, 
                                     ch_location, mg_lorentz, mg_couplings)

def add_vertices_to_calchep_files(ch_location, new_vertices, param_dict, lorentz_dict, 
                                  part_dict, mg_parts, mg_lorentz, mg_couplings):
    """
    Creates copies of the original CalcHEP files, and adds any new vertices found in the
    MadGraph files to them.
    """

    new_ch_loc = os.path.abspath(ch_location) + "_ufo2mdl"
    if os.path.isdir(new_ch_loc):
        print("\n{} already exists - deleting...").format(new_ch_loc)
        shutil.rmtree(new_ch_loc)        

    shutil.copytree(ch_location, new_ch_loc)

    # Firstly do the vertices
    for i in range(len(new_vertices)):
        vertex = new_vertices[i]
        write_ch_vertex(vertex, param_dict, lorentz_dict, part_dict, 
                        mg_parts, ch_location, mg_lorentz, mg_couplings)

    # Add Maux to the vars file, if we need auxiliary particles...
    with open(new_ch_loc+"/vars1.mdl", 'a') as f:
        f.write("Maux           |1            |Auxiliary mass parameter                                                  | \n")
        f.write("half           |0.5          |One divided by two                                                        | \n")

    ## SB 8/1/20 REMOVED AS OF NOW
    ## Go through the optimisations and add the new definitions
    #global optimisations
    # with open(new_ch_loc+"/func1.mdl", 'a') as f:
    #     for k, v in sorted(iteritems(optimisations)):
    #         towrite = "{0: <15}".format(k) + "|" + v + "  % Optimisation by GUM\n"
    #         f.write(towrite)

    print("New CalcHEP vertices added to {}.").format(new_ch_loc)

"""
Main functions
"""

def convert_mg_to_ch(mg_location):
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

    param_dict = {}
    lorentz_dict = {} 
    part_dict = {}

    print("Done.")

    print("Writing CalcHEP output...")

    write_calchep_files(params, verts, param_dict, 
                        lorentz_dict, part_dict, parts)

    print("Done.")
    
    
def compare_mg_and_ch(mg_location, ch_location):
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
        for x in mg_full:
            print(x.name)
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
    # different. They probably won't be if they've been generated by the same tool.
    part_dict = get_mg_ch_parts_dict(mg_parts, ch_parts)

    print("All particles are consistent...")
    print("Checking vertices...")

    lorentz_dict = get_mg_ch_lorentz_dict(mg_lorentz)

    all_present = compare_vertices(mg_verts, mg_parts, mg_params, mg_lorentz, 
                                   mg_couplings, ch_verts, ch_parts, ch_params,
                                   ch_location, param_dict, lorentz_dict, part_dict)

    print("Done!")

    # If they're all there then return the original CalcHEP location
    if all_present:
        print("All vertices consistent between CH and MG.")
        return ( os.path.abspath(ch_location) )
    else:
        print("Added new vertices to CalcHEP files.")
        return ( os.path.abspath(ch_location) + "_ufo2mdl" )
    return 
    
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
        convert_mg_to_ch(sys.argv[1])
        
    elif len(sys.argv) == 3:
        compare_mg_and_ch(sys.argv[1], sys.argv[2])
