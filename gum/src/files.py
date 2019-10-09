"""
Master module containing all routines for finding, creating,
amending, and reading files.
"""

import os
import re
import numpy as np
import yaml
from distutils.dir_util import remove_tree
from collections import defaultdict

from setup import *

def remove_tree_quietly(path):
    """
    Deletes a directory if it exists
    """
    if os.path.exists(path):
        remove_tree(path)

def mkdir_if_absent(path):
    """
    Makes a new directory if the proposed path doesn't yet exist
    """
    try:
        os.makedirs(path)
    except OSError:
        if not os.path.isdir(path):
            raise

def full_filename(filename, module):
    """
    Formats a gambit file correctly based on the filename, the
    module, and whether it is a header file.
    """

    # strip leading & trailing slashes
    filename.strip('/')
    module.strip('/')

    path = ""
    if filename.endswith(".hpp"):
        path = "include/gambit/{0}/".format(module)
    elif filename.endswith(".cpp"):
        path = "src/"
    elif filename.endswith(".py"):
        path = "scripts/"
    elif filename.endswith(".dif"):
        path = "patches/"
    else:
        path = ""

    location = "../{0}/{1}{2}".format(module, path, filename)
    return location

def find_file(filename, module):
    """
    Tries to find a file in a specified module.
    """

    location = full_filename(filename, module)

    if os.path.exists(location):
        return True
    else:
        return False

def find_capability(capability, module, filename=None):
    """
    Tries to find a CAPABILITY in a specified rollcall file.
    """

    module.strip('/')

    if not filename:
        filename = "{0}_rollcall.hpp".format(module)
    else:
        filename.strip('/')

    location = full_filename(filename, module)
    if find_file(filename, module):
        pass
    else:
        raise GumError(("\n\nCannot find capability {0} in rollcall header file"
                        " {1} as the file does not exist!!").format(capability,
                                                                    location))

    lookup = "#define CAPABILITY " + capability

    with open(location) as f:
        for num, line in enumerate(f, 1):
            if lookup in line:
                return True, num

    return False, 0

def amend_rollcall(capability, module, contents, reset_dict, filename=None):
    """
    Adds a new FUNCTION to an existing CAPABILITY in a rollcall header.
    """

    module.strip('/')

    if not filename:
        filename = "{0}_rollcall.hpp".format(module)
    else:
        filename.strip('/')

    location = full_filename(filename, module)
    if find_file(filename, module):
        pass
    else:
        raise GumError(("\n\nCannot find capability {0} in rollcall header file"
                        " {1} as the file does not exist!!").format(capability,
                                                                    location))

    lookup = "#define CAPABILITY " + capability

    found = False

    with open(location) as f:
        for num, line in enumerate(f, 1):
            if lookup in line:
                found = True
                break

    lookup = "#undef CAPABILITY"

    # Found the capability -> find the #undef
    if found == True:
        with open(location) as f:
            for i in xrange(num):
                f.next()
            for no, line in enumerate(f, 1+num):
                if lookup in line:
                    break
            amend_file(filename, module, contents, no-1, reset_dict)
    else:
        raise GumError(("\n\nCapability {0} not found in "
                        "{1}!").format(capability, filename))

def find_function(function, capability, module, filename=None):
    """
    Tries to find a FUNCTION in a specified rollcall header file.
    """

    module.strip('/')

    if not filename:
        filename = "{0}_rollcall.hpp".format(module)
    else:
        filename.strip('/')

    location = full_filename(filename, module)

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

def find_string(filename, module, string):
    """
    Tries to find a generic string in a given file.
    """

    location = full_filename(filename, module)

    if find_file(filename, module):
        pass
    else:
        raise GumError(("\n\nCannot find string {0} in the file"
                        " {1} as the file does not exist!!").format(string,
                                                                    location))

    with open(location) as f:
        for num, line in enumerate(f, 1):
            if string in line:
                return True, num

    return False, 0

def write_file(filename, module, contents, reset_dict):
    """
    Writes a file in a specified location.
    """

    location = full_filename(filename, module)
    location_parts = os.path.split(location)

    if find_file(filename, module):
        raise GumError(("\n\nTried to write file " + location +
                        ", but it already exists."))

    reset_dict['new_files']['files'].append(location)

    # Save a temp mug file, incase something goes wrong.
    drop_mug_file("mug_files/temp.mug", reset_dict)

    # Make the directory if it doesn't exist
    mkdir_if_absent(location_parts[0])
    # Create new file
    open(location, 'w').write(contents)

    print("File {} successfully created.".format(location))

def copy_file(filename, module, output_dir, reset_dict, existing=True):
    """
    Copies an output file in a specified location.
    """
    import shutil

    location = full_filename(filename, module)
    location_parts = os.path.split(location)
    GUM_version = output_dir + '/' + location[3:]

    # Todo - do something with this. It's not doing anything right now.
    if existing:
        reset_dict['copied_amended_files']['files'].append(location)
    else:
        reset_dict['new_files']['files'].append(location)

    # Save a temp mug file, incase something goes wrong.
    drop_mug_file("mug_files/temp.mug", reset_dict)

    mkdir_if_absent(location_parts[0])
    # Copy the file
    shutil.copy(GUM_version, location)

    print("Copied "+GUM_version+" to "+location+".")

def delete_file(filename, module):
    """
    Deletes a file in a specified location.
    """

    location = full_filename(filename, module)

    if find_file(filename, module):
        os.remove(location)
        print("File {} successfully removed.".format(location))

def amend_file(filename, module, contents, line_number, reset_dict):
    """
    Amends a file in a specified location with 'contents', starting
    from a given line number.
    """

    location = full_filename(filename, module)

    # Catch an error code
    if line_number == -1:
        raise GumError(("Error in amend_file routine. Received line_number"
                        " of " + line_number + " to write to the file "
                        "" + location + ". I think something is wrong."))

    if not find_file(filename, module):
        raise GumError(("\n\nERROR: Tried to amend file " + location +
                        ", but it does not exist."))

    # Check there's not already an identical entry - happens sometimes!
    present = False
    if location in reset_dict['amended_files']:
        if contents in reset_dict['amended_files'][location]:
            present = True

    if not present:
        reset_dict['amended_files'][location].append(contents)

    # Save a temp mug file, incase something goes wrong.
    drop_mug_file("mug_files/temp.mug", reset_dict)

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

def add_capability(module, capability, function, reset_dict,
                   returntype, filename=None, dependencies=None, 
                   allowed_models=None, backend_reqs=None):
    """
    Finds a capability in a given module. If it already exists, then add 
    a new function entry to it. If it is not found, then add the 
    capability too.
    """

    cap_exists, cap_line = find_capability(capability, module, filename)

    # Get the function signature for the rollcall
    func = write_function(function, returntype, dependencies, 
                          allowed_models, backend_reqs)

    # If the capability is there...
    if cap_exists:
        # Check the function isn't already there too
        func_exists, func_line = find_function(function, capability, 
                                               module, filename)

        if func_exists:
            raise GumError(("The function "  + function + " already exists "
                            "in the capability " + capability + " in the "
                            "module " + module + "."))

        # If there's already the capability, then the only thing we need 
        # is the function. Add it to the existing capability.
        amend_file(filename, module, func, cap_line+1, reset_dict)
        # The +1 takes into account the START_CAPABILITY line. Hopefully 
        # being a bit lazy here doesn't come back to bite me, but I doubt it.

    # Otherwise - add some capability tags either side and write it to the 
    # end of the file.
    else:
        contents = (
            "  #define CAPBILITY {0}\n"
            "  START_CAPABILITY\n"
            "  \n"
            "{1}"
            "  \n"
            "  #undef CAPABILITY\n"
        ).format(capability, func)

        n = -1
        location = full_filename(filename, module)
        with open(location) as f:
            for num, line in enumerate(f, 1):
                if "#undef CAPABILITY" in line: n = num
        amend_file(filename, module, contents, n+1, reset_dict)


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
        # List
        if isinstance(allowed_models, list):
            extras += "ALLOW_MODELS({0})\n".format(', '.join(allowed_models))
        else:
            extras += "ALLOW_MODELS({0})\n".format(allowed_models)

    if backend_reqs:
        for i in np.arange(len(backend_reqs)):
            bereq = backend_reqs[i]
            if not isinstance(bereq, BackendReq):
                raise GumError(("\n\nBackend requirement at position " + i +
                                "not passed as instance of class BackendReq."))

            add = ", ".join([bereq.capability,
                             "({0})".format(', '.join(bereq.tags)),
                             bereq.cpptype])

            if not bereq.var:
                add += ", ({0})".format(', '.join(bereq.args))

            extras += "BACKEND_REQ({0})\n".format(add)

    towrite = (
            "#define FUNCTION {0}\n"
            "  START_FUNCTION({1})\n"
            "{2}"
            "#undef FUNCTION\n"
    ).format(function, returntype, dumb_indent(2, extras))

    return dumb_indent(4, towrite)

def add_new_model_to_function(filename, module, capability, function, 
                              model_name, reset_dict, pattern="ALLOW_MODELS"):
    """
    Adds a new entry to the ALLOW_MODELS macro for a given (pre-existing)
    CAPABILITY and FUNCTION. Pattern can be overwritten by something else
    to match e.g. ALLOW_MODEL_DEPENDENCES
    """

    location = full_filename(filename, module)
    temp_location = location + "_temp"

    # Check the capability + function exist
    exists, num = find_function(function, capability, module, filename)
    if not exists:
        raise GumError(("Could not find function {0} in capability "
                        "{1} in {2}").format(function, capability, location))

    counter = 0
    take_it_slow = False
    modellist = ""
    adding_to_modellist = False
    done = False

    with open(location, 'r') as f, open(temp_location, 'w+') as g:
        # Write everything up to the function
        for line in f:
            counter += 1
            if counter > num and not done:
                # Be cool if we are nearby
                take_it_slow = True
            # If we're not nearby just go for it
            if not take_it_slow:
                g.write(line)
            # If we're nearby then go through line-by-line
            else:
                if pattern in line:
                    # Start adding to the list of models
                    adding_to_modellist = True
                elif not adding_to_modellist: 
                    g.write(line)
                if adding_to_modellist and not done:
                    modellist += line
                    # End of macro
                    if ")" in line: 
                        # Check there's no more than ten here. 
                        numentries = len(modellist.split(','))
                        # Add the model name to the end if there's
                        # less than ten (macro maximum)
                        if numentries < 10:
                            # Only replacing first occurence of a ")"
                            modellist = re.sub(r'\)', ", "+model_name+')', 
                                        modellist, 1) 
                        # Otherwise add a new entry altogether
                        else:
                            newentry = ""
                            # If it's something like MODEL_GROUP(group2(...))
                            # make a new group3 and add an 
                            # ALLOW_MODEL_COMBINATION(group1, group3)
                            # TODO check to see if group(n+1) already exists.
                            m = re.search(r'MODEL_GROUP\(group(\d)', pattern)
                            if m:
                                p = m.group(1) # Index
                                newentry = (
                                    "\n      MODEL_GROUP(group{0}({1}))"
                                    "\n      ALLOW_MODEL_COMBINATION(group1,group{0})"
                                ).format((int(p)+1), model_name)
                            # Write the new entry *before* the old one. This
                            # way it will be appended to first by GUM.
                            else:
                                newentry = "\n      {0}({1})".format(pattern, model_name)
                            g.write(newentry)
                        g.write(modellist)
                        adding_to_modellist = False
                        take_it_slow = False
                        done = True

    # Entry for the mug file to parse
    entry = location+'|'+capability+'|'+function+'|'+pattern
    # Add to the reset dictionary
    reset_dict['new_models'][entry].append(model_name)

    os.remove(location)
    os.rename(temp_location, location)

    print("Model {} added to capability {}.".format(model_name, capability))


def blame_gum(message):
    """
    Writes function to dump at the beginning of a new GAMBIT file.
    Blames GUM. Takes a message to describe the new file.
    """

    if not message.startswith("///"):
        message = "/// " + message

    towrite = (
            "//   GAMBIT: Global and Modular BSM Inference Tool\n"
            "//   *********************************************\n"
            "///  \\file\n"
            "///\n"
            + message +
            "\n"
            "///\n"
            "///  Authors (add name and date if you modify):    \n"
            "///       *** Automatically created by GUM ***     \n"
            "///                                                \n"
            "///  \\author The GAMBIT Collaboration             \n"
            "///  \date " +
            datetime.datetime.now().strftime("%I:%M%p on %B %d, %Y") +
            "\n"
            "///                                                \n"
            "///  ********************************************* \n"
            "\n"
    )

    return towrite

def dumb_indent(numspaces, text):
    """
    Indents input text by specified number of spaces.
    """

    lines = text.splitlines(True)
    out = ""
    for line in lines:
       for i in range(0, numspaces):
           out += " "
       out += line

    return out

def indent(text, numspaces=0):
    """
    Increases indents of text every time it detects an open
    curly brace {, and decreases every time it detects a closed
    curly brace }. Does nothing if there is a {}.
    """

    lines = text.splitlines(True)
    out = ""

    for line in lines:


        if numspaces < 0:
            raise GumError(("\n\nTried to indent a negative number of spaces."
                            "Please check for rogue braces."))

        num_open = line.count("{")
        num_closed = line.count("}")
        if num_open != num_closed:
            numspaces -= (num_closed*2)
        for i in range(0, numspaces):
            out += " "
        if num_open != num_closed:
            numspaces += (num_open*2)
        out += line

    return out

def revert(reset_file):
    """
    Go back to the previous save point.
    """

    print("GUM called in reset mode.")

    if not reset_file.endswith('.mug'):
        raise GumError(("\n\n\tPlease use a generated .mug "
                       "file if using the reset function."))

    with open(reset_file, "r") as f:

        try:
            data = yaml.safe_load(f)
        except yaml.YAMLError as exc:
            print(exc)
            return

        # The files GUM wrote as new.
        # GUM can just simply delete these.
        new_files = data['new_files']['files']


        for i in new_files:
            if os.path.isfile(i):
                print("Deleting {0}...".format(i))
                os.remove(i)
            else:
                print("Tried deleting {0}, but it seems to have already been removed.".format(i))

        # The files that existed previously, that GUM added stuff to.
        # These are a little more annoying to deal with.
        amended_files = data['amended_files']

        # We want to match *strings* and not line numbers or anything like that.
        # This way, there is no order needed to perform resets in.

        for filename, v in amended_files.iteritems():

            print("Amending {0}...".format(filename))

            temp_file = filename + "_temp"
            with open(filename, 'r') as original_file:
                text = original_file.read()
                lines = text.splitlines()

            # This is the YAML node for the strings needing to be deleted
            for string in v:

                # If we find it, kill it
                if string in text:
                    text = text.replace(string, '')
                # Otherwise assume someone has messed with it
                else:
                    print(("Tried deleting the following entry from "
                           "{0}, but I could not find it:\n\n{1}\n\n"
                           "Perhaps it has already been removed."
                           ).format(filename, string))

            # Write the new amended file
            new_file = open(temp_file, 'w')
            new_file.write(text)

            os.remove(filename)
            os.rename(temp_file, filename)

        # Now go through those amendments that are adding new models to an existing
        # ALLOW_MODELS (or similar) macro for a given capability
        amended_capabilities = data['new_models']

        for loc_cap_func_pattern, model in amended_capabilities.iteritems():
            
            location, capability, function, pattern = loc_cap_func_pattern.split('|')
            module = location.split('/')[1]

            print(("Removing model from capability {0}; function {1}; in {2}...")
                   .format(capability, function, module))

            temp_file = location + "_temp"
            
            exists, num = find_function(function, capability, module)
            
            counter = 0
            done = False
            take_it_slow = False
            modellist = ""
            taking_from_modellist = False
            
            with open(location, 'r') as f, open(temp_file, 'w+') as g:
                for line in f:
                    counter += 1
                    if counter > num and not done: 
                        take_it_slow = True
                    if not take_it_slow: 
                        g.write(line)
                    else:
                        if pattern in line: 
                            taking_from_modellist = True 
                        elif not taking_from_modellist: g.write(line)
                        if taking_from_modellist and not done:
                            modellist += line
                            if ")" in line: 
                                # Take the model out of the list
                                modellist = re.sub(', '+model[0], '', modellist)
                                g.write(modellist)
                                taking_from_modellist = False
                                take_it_slow = False
                                done = True

            os.remove(location)
            os.rename(temp_file, location)

    return

def check_for_existing_entries(model_name, darkbit, colliderbit, output_opts):
    """
    Checks for existing entries within GAMBIT, for a new model.
    """

    # Models entries
    if ( find_file("models/" + model_name + ".hpp", "Models") or
         find_file("SpectrumContents/" + model_name + ".cpp", "Models") or 
         find_file("SimpleSpectra/" + model_name + "SimpleSpec" + ".hpp", "Models") or 
         find_file("SpecBit_" + model_name + ".cpp", "SpecBit") or 
         find_file("SpecBit_" + model_name + "_rollcall.hpp", "SpecBit")
       ):
        raise GumError(("Model {0} already exists in the Model Hierarchy.").format(model_name))

    if darkbit:
        if find_file(model_name + ".cpp", "DarkBit"):
            raise GumError(("Model {0} already exists in DarkBit.").format(model_name))
        if output_opts.mo:
            ver = "3.6.9.2"
            f = "frontends/MicrOmegas_{0}_{1}".format(model_name, ver.replace('.','_'))
            if ( find_file(f+".cpp", "Backends") or 
                 find_file(f+".hpp", "Backends")
               ):
                raise GumError(("MicrOmegas entry already exists for model {0}").format(model_name))

def drop_mug_file(mug_file, contents):
    """
    Drops a .mug (reset) file from the reset contents saved by GUM.
    """

    d = dict(contents)

    if 'new_files' in d:
        new_files = dict(d['new_files'])
    else:
        new_files = {}
    if 'amended_files' in d:
        amended_files = dict(d['amended_files'])
    else:
        amended_files = {}
    if 'new_models' in d:
        new_models = dict(d['new_models'])
    else:
        new_models = {}
        
    new_contents = {'new_files': new_files, 'amended_files': amended_files, 
                    'new_models' : new_models}

    with open(mug_file, 'w') as f:
        yaml.dump(new_contents, f, default_flow_style=False)


def drop_yaml_file(model_name, model_parameters, add_higgs, reset_contents):
    """
    Drops an example YAML file with all decays of a new model
    added.
    """

    towrite = (
        "##########################################################################\n"
        "## GAMBIT configuration for a random scan of the {0} model\n"
        "##\n"
        "## Includes simply the decays of the new model.\n"
        "##########################################################################\n"
        "\n"
        "\n"
        "Parameters:\n"
        "\n"
        "  # SM parameters.\n"
        "  StandardModel_SLHA2: !import include/StandardModel_SLHA2_defaults.yaml\n"
        "\n"
    ).format(model_name)

    if add_higgs:
        towrite+= (
            "  StandardModel_Higgs:\n"
            "    mH: 125.09\n"
            "\n"
        )

    towrite += ("  {0}:\n").format(model_name)

    # Don't want the SM-like Higgs mass a fundamental parameter
    bsm_params = [x for x in model_parameters if x.name != 'h0_1' and x.sm == False]
    params = []

    for i in bsm_params:
        if i.shape == 'scalar' or i.shape == None: params.append(i.gb_in)
        elif re.match("m[2-9]x[2-9]", i.shape): 
            size = int(i.shape[-1])
            for j in xrange(size):
                for k in xrange(size):
                    params.append(i.gb_in + str(j+1) + 'x' + str(k+1))
        elif re.match("v[2-9]", i.shape): 
            size = int(i.shape[-1])
            for j in xrange(size):
                params.append(i.gb_in + str(j+1))

    # No double counting (also want to preserve the order)
    norepeats = []
    [norepeats.append(i) for i in params if not i in norepeats]

    for i in norepeats:
        towrite += ("    {0}: 0.1\n").format(i)

    towrite += (
        "Priors:\n"
        "\n"
        "  # All the priors are simple for this scan, so they "
        "are specified directly in the Parameters section.\n"
        "\n"
        "\n"
        "Printer:\n"
        "\n"
        "  printer: hdf5\n"
        "\n"
        "  options:\n"
        "    output_file: \"{0}.hdf5\"\n"
        "    group: \"/{0}\"\n"
        "\n"
        "\n"
        "Scanner:\n"
        "\n"
        "  use_scanner: random\n"
        "\n"
        "  scanners:\n"
        "\n"
        "    random:\n"
        "      plugin: random\n"
        "      point_number: 10\n"
        "      like:  LogLike\n"
        "\n"
        "ObsLikes:\n"
        "\n"
        "  - purpose:      Observable\n"
        "    capability:   decay_rates\n"
        "    type:         DecayTable\n"
        "    printme:      false\n"
        "\n"
        "Rules:\n"
        "\n"
        "  # Choose to get decays from DecayBit proper, not from an SLHA file.\n"
        "  - capability: decay_rates\n"
        "    function: all_decays\n"
        "Logger:\n"
        "\n"
        "  redirection:\n"
        "    [Debug] : \"debug.log\"\n"
        "    [Default] : \"default.log\"\n"
        "    [DecayBit] : \"DecayBit.log\"\n"
        "    [PrecisionBit] : \"PrecisionBit.log\"\n"
        "    [ColliderBit] : \"ColliderBit.log\"\n"
        "    [SpecBit] : \"SpecBit.log\"\n"
        "    [Dependency Resolver] : \"dep_resolver.log\"\n"
        "\n"
        "KeyValues:\n"
        "\n"
        "  dependency_resolution:\n"
        "    prefer_model_specific_functions: true\n"
        "\n"
        "  likelihood:\n"
        "    model_invalid_for_lnlike_below: -5e5\n"
        "    model_invalid_for_lnlike_below_alt: -1e5\n"
        "\n"
        "  default_output_path: \"runs/{0}/\"\n"
        "\n"
        "  debug: false\n"
        "\n"
    ).format(model_name)

    write_file(model_name + '.yaml', 'yaml_files', towrite, reset_contents)

def write_config_file(outputs, model_name, reset_contents):
    """
    Drops a configuration file, which will build the correct backends, 
    and then GAMBIT, in the correct order.
    """

    towrite = (
        "cd ../build\n"
        "cmake ..\n"
        "make -j4"
    )

    if outputs.pythia:
        towrite += " pythia_{0}\n".format(model_name.lower())

    if outputs.mo:
        towrite += " micromegas_{0}\n".format(model_name)

    if outputs.spheno:
        towrite += " spheno_{0}".format(model_name.lower())

    if outputs.vev:
        towrite += " vevacious"

    if outputs.ch:
        towrite += " calchep"

    # TODO : flexiblesusy.     
    towrite += (
        "\n"
        "cmake ..\n"      # Have to cmake here because of Pythia headers.
        "make -j4 gambit\n"
    )

    write_file(model_name + '_config.sh', 'gum', towrite, reset_contents)