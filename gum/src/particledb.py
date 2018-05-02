"""
Contains scripts for manipulating the particle database. This includes
scraping the information about existing particles as well as adding
new particles to the database.
"""

import yaml
import os
import itertools
import sys


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

        # Filter out comments, etc. Only care about those with a PDG_context pair and a name.
        if 'PDG_context' in particles[i]:

            # Add to GAMBIT dict.
            gambit_pdg_codes[particles[i]['name']] = \
                particles[i]['PDG_context'][0]

            # Add the conjugate as well.
            if 'conjugate' in particles[i]:
                gambit_pdg_codes[particles[i]['conjugate']] = \
                - particles[i]['PDG_context'][0]

            # Is there a DecayBit entry?
            if 'DecayBit' in particles[i]:

                # If particle has decays
                if particles[i]['DecayBit']['Decays'] == True:

                    # If particle has DecayBit name specified
                    if 'name' in particles[i]['DecayBit']:
                        decaybit_dict[particles[i]['DecayBit']['name']] = \
                        particles[i]['PDG_context'][0]
                    else:
                        decaybit_dict[particles[i]['name']] = \
                        particles[i]['PDG_context'][0]

                    # Likewise for conjugates
                    if 'conjugate' in particles[i]:
                        if 'conjugate' in particles[i]['DecayBit']:
                            decaybit_dict[
                                particles[i]['DecayBit']['conjugate']] = \
                            - particles[i]['PDG_context'][0]
                        else:
                            decaybit_dict[particles[i]['conjugate']] = \
                            - particles[i]['PDG_context'][0]

    # Now the particle sets
    sets = data['StandardModel']['Sets'] + data['OtherModels']['Sets']

    for i in range(0, len(sets)):

        if 'PDG_context' in sets[i]:
            for j in range(0, len(sets[i]['PDG_context'])):
                gambit_pdg_codes[sets[i]['name'] + "_" + str(j + 1)] = \
                sets[i]['PDG_context'][j][0]

            if 'conjugate' in sets[i]:
                for j in range(0, len(sets[i]['PDG_context'])):
                    gambit_pdg_codes[
                        sets[i]['conjugate'] + "_" + str(j + 1)] = \
                    - sets[i]['PDG_context'][j][0]

            # If a particle set has Decays
            if 'DecayBit' in sets[i]:

                # Do we care about the Decays?
                if sets[i]['DecayBit']['Decays'] == True:

                    # If the name is specified, use these
                    if 'name' in sets[i]['DecayBit']:

                        # If each particle in the set has a different DecayBit
                        # rollcall name (e.g. Higgs, h0_2) then parse on a
                        # case-by-case
                        if type(sets[i]['DecayBit']['name']) is list:
                            for j in range(0, len(sets[i]['PDG_context'])):
                                decaybit_dict[sets[i]['DecayBit']['name'][j]] = \
                                sets[i]['PDG_context'][j][0]

                        # Otherwise iterate over the set (particle_1,
                        # particle_2, etc.)
                        elif type(sets[i]['DecayBit']['name']) is str:
                            for j in range(0, len(sets[i]['PDG_context'])):
                                decaybit_dict[
                                    sets[i]['DecayBit']['name'] + "_" + str(
                                        j + 1)] = sets[i]['PDG_context'][j][0]

                    # If the name is unspecified, use the default
                    else:
                        if type(sets[i]['name']) is list:
                            for j in range(0, len(sets[i]['PDG_context'])):
                                decaybit_dict[sets[i]['name'][j]] = \
                                sets[i]['PDG_context'][j][0]

                        elif type(sets[i]['name']) is str:
                            for j in range(0, len(sets[i]['PDG_context'])):
                                decaybit_dict[
                                    sets[i]['name'] + "_" + str(j + 1)] = \
                                sets[i]['PDG_context'][j][0]

                    # Conjugates
                    if 'conjugate' in sets[i]:

                        # If specified names, use these
                        if 'conjugate' in sets[i]['DecayBit']:

                            if type(sets[i]['DecayBit']['conjugate']) is list:
                                for j in range(0, len(sets[i]['PDG_context'])):
                                    decaybit_dict[
                                        sets[i]['DecayBit']['conjugate'][j]] = \
                                    - sets[i]['PDG_context'][j][0]
                            elif type(sets[i]['DecayBit']['conjugate']) is str:
                                for j in range(0, len(sets[i]['PDG_context'])):
                                    decaybit_dict[sets[i]['DecayBit'][
                                                      'conjugate'] + "_" + str(
                                        j + 1)] = -sets[i]['PDG_context'][j][0]

                        # Otherwise default
                        else:

                            if type(sets[i]['conjugate']) is list:
                                for j in range(0, len(sets[i]['PDG_context'])):
                                    decaybit_dict[sets[i]['conjugate'][j]] = - \
                                    - sets[i]['PDG_context'][j][0]
                            elif type(sets[i]['conjugate']) is str:
                                for j in range(0, len(sets[i]['PDG_context'])):
                                    decaybit_dict[
                                        sets[i]['conjugate'] + "_" + str(
                                            j + 1)] = - \
                                    sets[i]['PDG_context'][j][0]

    return gambit_pdg_codes, decaybit_dict

def check_all_particles_present(partlist, gambit_pdg_codes):
    """
    Checks all particles exist in the particle_database.yaml.
    """
    

def add_new_particleDB_entry(particle):
    """
    Adds a new particle to the particle_database.yaml file.
    """
