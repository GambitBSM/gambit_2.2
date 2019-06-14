//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  Standalone executable to test basic Spectrum
///  and SpectrumContents classes 
///
///  Not built with GAMBIT cmake system. 
///  Compile like this:
///
///  > g++ -Wfatal-errors --std=c++14 -o test_spectrum_objects \
-fopenmp \
-Lcontrib/yaml-cpp-0.6.2 \
-lyaml-cpp \
-IModels/include \
-IElements/include \
-IUtils/include \
-ILogs/include \
-Icontrib/mkpath/include \
-Icontrib/yaml-cpp-0.6.2/include \
-Icontrib/slhaea/include \
-Icmake/include \
contrib/mkpath/src/mkpath.c \
Elements/src/spectrum.cpp \
Elements/src/slhaea_helpers.cpp \
Utils/src/exceptions.cpp \
Utils/src/standalone_error_handlers.cpp \
Utils/src/util_functions.cpp \
Logs/src/logger.cpp \
Logs/src/logmaster.cpp \
Logs/src/logging.cpp \
Models/src/spectrum_contents.cpp \
Models/src/particle_database.cpp \
Models/src/partmap.cpp \
Models/src/spectrum_contents.cpp \
Models/src/SpectrumContents/DiracSingletDM.cpp \
Models/src/SpectrumContents/MajoranaSingletDM.cpp \
Models/src/SpectrumContents/MDM.cpp \
Models/src/SpectrumContents/MSSM.cpp \
Models/src/SpectrumContents/ScalarSingletDM.cpp \
Models/src/SpectrumContents/SM.cpp \
Models/src/SpectrumContents/SMHiggs.cpp \
Models/src/SpectrumContents/SM_slha.cpp \
Models/src/SpectrumContents/VectorSingletDM.cpp \
SpecBit/tests/spectrum_object_tests.cpp 


//
//    -IModels/include
//    -IElements/include
//    -IUtils/include
///
///  *********************************************
///
///  Authors (add name and date if you modify):
///
///  \author Ben Farmer
///          (benjamin.farmer@imperial.ac.uk)
///  \date 2019 June
///
///  *********************************************

#include "gambit/Utils/standalone_error_handlers.hpp"
#include "gambit/Utils/standalone_utils.hpp"
#include "gambit/Utils/static_members.hpp"
#include "gambit/Logs/logger.hpp"
#include "gambit/Models/SpectrumContents/spectrum_contents.hpp"
#include "gambit/Models/SpectrumContents/RegisteredSpectra.hpp"
#include "gambit/Elements/spectrum.hpp"

using namespace Gambit;

int main(int argc, char* argv[])
{
    std::cout<<"Creating MSSM Spectrum Contents object"<<std::endl; 
    // Create MSSM spectrum contents object
    SpectrumContents::MSSM mssm; 

    std::cout<<"Writing template SLHA file corresponding to required spectrum contents: mssm_template.slha"<<std::endl;
    // Create template SLHA file complying with MSSM contents definition
    // (Won't include stuff the Spectrum objects don't care about, like MODSEL blocks and whatnot)
    mssm.create_template_SLHA_file("mssm_template.slha");
}
