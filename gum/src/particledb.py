"""
Contains scripts for manipulating the particle database. This includes
scraping the information about existing particles as well as adding
new particles to the database.
"""

import yaml
import os
import itertools
import sys

from setup import *

def get_gambit_particle_pdg_dict():
    """
    Return dict of GAMBIT particle names & PDG-context pairs from
    particle_database.yaml.
    """

    # TODO: Error checking of lengths of arrays etc.

    with open("./../config/particle_database.yaml", "r") as f:

        try:
            data = yaml.load(f)
        except yaml.YAMLerror as exc:
            print(exc)

    gambit_pdg_codes = {}
    decaybit_dict = {}

    # Firstly deal with particles
    particles = data['StandardModel']['Particles'] \
              + data['OtherModels']['Particles']

    for i in range(0, len(particles)):
    
        part = particles[i]

        # Filter out comments, etc. Only care about those with a 
        # PDG_context pair (plus a name).
        if 'PDG_context' in part:
            
            # Only add mass eigenstates
            if not (part['PDG_context'][1] == 0):
                continue
            
            PDG_code = part['PDG_context'][0]
        
            # If the particle is already in the dictionary, don't add it again
            if PDG_code in gambit_pdg_codes.values():
                continue
                
            # Add to GAMBIT dict.
            gambit_pdg_codes[part['name']] = PDG_code

            # Add the conjugate as well.
            if 'conjugate' in part:
                gambit_pdg_codes[part['conjugate']] = -PDG_code

            # Is there a DecayBit entry?
            if 'DecayBit' in part:

                # If particle has decays
                if part['DecayBit']['Decays'] == True:

                    # If particle has DecayBit name specified
                    if 'name' in part['DecayBit']:
                        decaybit_dict[part['DecayBit']['name']] = PDG_code
                    else:
                        decaybit_dict[part['name']] = PDG_code

                    # Likewise for conjugates
                    if 'conjugate' in part:
                        if 'conjugate' in part['DecayBit']:
                            decaybit_dict[part['DecayBit']['conjugate']] = -PDG_code
                        else:
                            decaybit_dict[part['conjugate']] = -PDG_code

    # Now the particle sets
    sets = data['StandardModel']['Sets'] + data['OtherModels']['Sets']

    for i in range(0, len(sets)):
    
        cset = sets[i]

        if 'PDG_context' in cset:
        
            PDG_codes = [x[0] for x in cset['PDG_context']]
            
            # Only add mass eigenstates.
            # TODO - figure out what to do with e.g. nu_RPV.
            if not (cset['PDG_context'][0][1] == 0):
                continue
            
            for PDG_code in PDG_codes:
                        
                # Don't add duplicates. Shouldn't be any, anyway...
                if PDG_code in gambit_pdg_codes.values():
                    continue
                    
                name = cset['name']
        
                for j in range(0, len(cset['PDG_context'])):
                    gambit_pdg_codes[name + "_" + str(j+1)] = cset['PDG_context'][j][0]
    
                if 'conjugate' in cset:
                    for j in range(0, len(cset['PDG_context'])):
                        gambit_pdg_codes[cset['conjugate'] + "_" + str(j+1)] = -cset['PDG_context'][j][0]
    
                # If a particle set has Decays
                if 'DecayBit' in cset:
    
                    # Do we care about the Decays?
                    if cset['DecayBit']['Decays'] == True:
    
                        # If the name is specified, use these
                        if 'name' in cset['DecayBit']:
    
                            # If each particle in the set has a different DecayBit
                            # rollcall name (e.g. Higgs, h0_2) then parse on a
                            # case-by-case
                            if type(cset['DecayBit']['name']) is list:
                                for j in range(0, len(cset['PDG_context'])):
                                    decaybit_dict[cset['DecayBit']['name'][j]] = cset['PDG_context'][j][0]
    
                            # Otherwise iterate over the set (part_1, part_2, etc.)
                            elif type(cset['DecayBit']['name']) is str:
                                for j in range(0, len(cset['PDG_context'])):
                                    decaybit_dict[cset['DecayBit']['name'] + "_" + str(j+1)] = cset['PDG_context'][j][0]
    
                        # If the name is unspecified, use the default
                        else:
                            if type(cset['name']) is list:
                                for j in range(0, len(cset['PDG_context'])):
                                    decaybit_dict[cset['name'][j]] = cset['PDG_context'][j][0]
    
                            elif type(cset['name']) is str:
                                for j in range(0, len(cset['PDG_context'])):
                                    decaybit_dict[cset['name'] + "_" + str(j+1)] = cset['PDG_context'][j][0]
    
                        # Conjugates
                        if 'conjugate' in cset:
    
                            # If specified names, use these
                            if 'conjugate' in cset['DecayBit']:
    
                                if type(cset['DecayBit']['conjugate']) is list:
                                    for j in range(0, len(cset['PDG_context'])):
                                        decaybit_dict[cset['DecayBit']['conjugate'][j]] = - cset['PDG_context'][j][0]
                                        
                                elif type(cset['DecayBit']['conjugate']) is str:
                                    for j in range(0, len(cset['PDG_context'])):
                                        decaybit_dict[cset['DecayBit']['conjugate'] + "_" + str(j+1)] = -cset['PDG_context'][j][0]
    
                            # Otherwise default
                            else:
    
                                if type(cset['conjugate']) is list:
                                    for j in range(0, len(cset['PDG_context'])):
                                        decaybit_dict[cset['conjugate'][j]] = - cset['PDG_context'][j][0]
                                elif type(cset['conjugate']) is str:
                                    for j in range(0, len(cset['PDG_context'])):
                                        decaybit_dict[cset['conjugate'] + "_" + str(j+1)] = -cset['PDG_context'][j][0]

    return gambit_pdg_codes, decaybit_dict

def check_all_particles_present(partlist, gambit_pdg_codes):
    """
    Checks all particles exist in the particle_database.yaml.
    """
    
    absent = []
    
    for i in xrange(len(partlist)):
        if not partlist[i].pdg() in gambit_pdg_codes.values():
            absent.append(partlist[i].pdg())
    
    if len(absent) == 0:
        print("All particles are in the GAMBIT database.")
        return True
    else:
        print(("\n\tThe following particles (by PDG code) are missing from the "
               "particle database: {0}. Please add them to "
               "../config/particle_database.yaml.\n").format(absent))
        return False
    

def add_new_particleDB_entry(particle):
    """
    Adds a new particle to the particle_database.yaml file.
    """
