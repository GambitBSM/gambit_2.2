#!/bin/bash

g++ -g -Wfatal-errors --std=c++14 -o test_spectrum_objects \
    -fopenmp \
    -IModels/include \
    -IElements/include \
    -IUtils/include \
    -ILogs/include \
    -I/usr/local/include \
    -Icontrib/slhaea/include \
    -Icmake/include \
    Elements/src/spectrum.cpp \
    Elements/src/slhaea_helpers.cpp \
    Utils/src/exceptions.cpp \
    Utils/src/standalone_error_handlers.cpp \
    Utils/src/util_functions.cpp \
    Utils/src/file_lock.cpp \
    Utils/src/version.cpp \
    Logs/src/logger.cpp \
    Logs/src/logmaster.cpp \
    Logs/src/logging.cpp \
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
    SpecBit/tests/spectrum_object_tests.cpp \
    -lyaml-cpp \
