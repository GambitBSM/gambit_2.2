"""
Contains all routines for Vevacious output from SARAH.
"""

from setup import *
from files import *

import shutil
import os

def copy_vevacious_files(model_name, vevdir):
    """
    Moves vevacious files from the SARAH output to the backend folder.
    """

    gb_target = "./../Backends/installed/vevacious/VevaciousPlusPlus/1.0/ModelFiles/" + model_name
    if not os.path.exists(gb_target):
        os.makedirs(gb_target)

    file = [f for f in os.listdir(vevdir) if f.endswith(".vin")]
    if len(file) > 1:
        raise GumError(("Too many .vin files in vevacious output."))
    if len(file) == 0:
        raise GumError(("No .vin files in vevacious output."))

    shutil.copyfile(vevdir + "/" + file[0], gb_target + "/" + model_name + ".vin")
    shutil.copyfile(vevdir + "/ScaleAndBlock.xml", gb_target + "/ScaleAndBlock.xml")

def write_vevacious_src(model_name, vevdir, spectrum, parameters):
    """
    Writes source code for a new vevacious model.
    To go in SpecBit/src/SpecBit_VS.cpp.
    """

    copy_vevacious_files(model_name, vevdir)

    towrite = (
        "/******************************/\n"
        "/* {0:^26} */\n"
        "/******************************/\n"
        "\n"
        "/// Tell GAMBIT which files to work with for the {0} model.\n"
        "void vevacious_file_location_{0}(map_str_str &result)\n"
        "{{\n"
        "namespace myPipe = Pipes::vevacious_file_location_{0};\n"
        "const Options& runOptions = *myPipe::runOptions;\n"
        "\n"
        "int rank;\n"
        "// Get mpi rank\n"
        "#ifdef WITH_MPI\n"
        "    MPI_Comm_rank(MPI_COMM_WORLD, &rank);\n"
        "#else\n"
        "    rank = 0;\n"
        "#endif\n"
        "\n"
        "// Creating string with rank number\n"
        "std::string rankstring = std::to_string(rank);\n"
        "\n"
        "// Getting the run folder for saving initialization files\n"
        "std::string inputspath = runOptions.getValue<std::string>(\"where_to_save_input\");\n"
        "result[\"inputspath\"] = inputspath;\n"
        "std::string modelfilesPath = inputspath + \"/ModelFiles/mpirank_\"+ rankstring + \"/\";\n"
        "\n"
        "// Get the path to the library\n"
        "std::string vevaciouslibpath = Backends::backendInfo().path_dir(\"vevacious\", \"1.0\");\n"
        "std::string vevaciouspath = vevaciouslibpath + \"/../\";\n"
        "\n"
        "result[\"ScaleAndBlockFileSource\"] = vevaciouspath + \"ModelFiles/\"{0}/ScaleAndBlock.xml\";\n"
        "result[\"ModelFileSource\"] = vevaciouspath + \"ModelFiles/{0}/{0}.vin\";\n"
        "result[\"ScaleAndBlockFile\"] = modelfilesPath + \"ScaleAndBlockFile.xml\";\n"
        "result[\"ModelFile\"] = modelfilesPath + \"ModelFile.vin\";\n"
        "}}\n"
        "\n"
        "// This function gives back the result for absolute stability, either \"Stable\""
        " or \"Metastable\".\n"
        "void check_stability_{0}(VevaciousResultContainer &result)\n"
        "{{\n"
        "namespace myPipe = Pipes::check_stability_{0};\n"
        "\n"
        "//static std::string inputspath =  *myPipe::Dep::make_vevaciousPlusPlus_inputs;\n"
        "static std::string inputspath =  *myPipe::Dep::init_vevacious;\n"
        "\n"
        "// Reset all member variables of VevaciousResultContainer to -1\n"
        "// to avoid that any value could be carried over from a previous calculated point\n"
        "result.reset_results();\n"
        "\n"
        "// Getting mpi rank\n"
        "int rank;\n"
        "#ifdef WITH_MPI\n"
        "            MPI_Comm_rank(MPI_COMM_WORLD, &rank);\n"
        "#else\n"
        "            rank = 0;\n"
        "#endif\n"
        "\n"
        "std::string rankstring = std::to_string(rank);\n"
        "\n"
        "std::string inputFilename = inputspath + \"/InitializationFiles/VevaciousPlusPlusObjectInitialization_mpirank_\"+ rankstring +\".xml\";\n"
        "vevacious_1_0::VevaciousPlusPlus::VevaciousPlusPlus vevaciousPlusPlus( inputFilename );\n"
        "\n"
        "// Get the spectrum object\n"
        "const Spectrum& fullspectrum = *myPipe::Dep::{2};\n"
        "const SubSpectrum& spectrumHE = fullspectrum.get_HE();\n"
        "\n"
        "// Get the SLHAea::Coll object from the spectrum.\n"
        "SLHAea::Coll slhaea = spectrumHE.getSLHAea(2);\n"
        "\n"
        "double scale = spectrumHE.GetScale();\n"
        "\n"
    ).format(model_name, vevdir, spectrum)

    # Now we go through each SLHA BLOCK that we have in the model, and use
    # the vevacious ReadLhaBlock routine from the SLHAea object, to pass the
    # potential parameters to vevacious.
    towrite +="// Here we start passing the parameters form the SLHAea::Coll object.\n"

    #towrite += write_readlhablock("YE", "m3x3")

    # print indent(towrite)
    #return indent(towrite)

    """
    TODO:

    - get all expected blocks from vevacious
    - need the following:
        - input parameters (SMINPUTS, Yukawas, etc.)
        - output parameters (GAUGE, mixings, some masses.)

    - function to return a list/dict of all BLOCKS and their entries. in and out. 
      needs the size of each block etc.

    """

    for param in parameters:

        print param.name, param.block, param.index

        # std::vector<std::pair<int,double>> gaugecouplings = 
        # { { 1 , SLHAea::to<double>(slhaea.at("GAUGE").at(1).at(1))  }, { 2, SLHAea::to<double>(slhaea.at("GAUGE").at(2).at(1)) }, { 3, SLHAea::to<double>(slhaea.at("GAUGE").at(3).at(1)) } };
        
        # vevaciousPlusPlus.ReadLhaBlock( "GAUGE", scale , gaugecouplings, 1 );

        # std::vector<std::pair<int,double>> Hmix = { { 1 , SLHAea::to<double>(slhaea.at("HMIX").at(1).at(1))},
        #                       { 101, SLHAea::to<double>(slhaea.at("HMIX").at(101).at(1))},
        #                       { 102, SLHAea::to<double>(slhaea.at("HMIX").at(102).at(1))},
        #                       { 103, SLHAea::to<double>(slhaea.at("HMIX").at(103).at(1))},
        #                       { 3, SLHAea::to<double>(slhaea.at("HMIX").at(3).at(1))}
        #                       };

        # vevaciousPlusPlus.ReadLhaBlock( "HMIX", scale , Hmix, 1 );
      
       
        # std::vector<std::pair<int,double>> minpar = {  
        #                       { 3, SLHAea::to<double>(slhaea.at("MINPAR").at(3).at(1))}
        #                       };
                            
        # vevaciousPlusPlus.ReadLhaBlock( "MINPAR", scale , minpar, 1 );

        #     std::vector<std::pair<int,double>> msu2 = { { 11 , SLHAea::to<double>(slhaea.at("MSU2").at(1,1).at(2))},
        #                                                 { 12, SLHAea::to<double>(slhaea.at("MSU2").at(1,2).at(2))},
        #                                                 { 13, SLHAea::to<double>(slhaea.at("MSU2").at(1,3).at(2))},
        #                                                 { 21, SLHAea::to<double>(slhaea.at("MSU2").at(2,1).at(2))},
        #                                                 { 22, SLHAea::to<double>(slhaea.at("MSU2").at(2,2).at(2))},
        #                                                 { 23, SLHAea::to<double>(slhaea.at("MSU2").at(2,3).at(2))},
        #                                                 { 31, SLHAea::to<double>(slhaea.at("MSU2").at(3,1).at(2))},
        #                                                 { 32, SLHAea::to<double>(slhaea.at("MSU2").at(3,2).at(2))},
        #                                                 { 33, SLHAea::to<double>(slhaea.at("MSU2").at(3,3).at(2))}
        #     };
                              
        # vevaciousPlusPlus.ReadLhaBlock( "MSU2", scale , msu2, 2 );


def write_readlhablock(block, shape, entries=[]):
    """
    TODO
    Writes the ReadLhaBlock routine for a given block.
    """

    towrite = ""
    pad = 40 + len(block)

    # Matrix case
    if shape.startswith('m'):

        # The size of the matrix
        size = int(shape.split('x')[-1])

        towrite += (
                "std::vector<std::pair<int,double>> {0} = {{ "
        ).format(block)

        # Now do each element of the matrix
        entries = []
        for i in range(size):
            for j in range(size):
                entry = (
                    "{{ {0}{1}, SLHAea::to<double>(slhaea.at(\"{2}\")"
                    ".at({0},{1}),at(2))}}"
                ).format(i+1, j+1, block)
                entries.append(entry)
        
        towrite += ",\n{:>{width}}".format("",width=pad).join(entries) +"\n};\n"
        towrite += "vevaciousPlusPlus.ReadLhaBlock(\"{0}\", scale, {0}, 2);\n".format(block)

    return towrite

