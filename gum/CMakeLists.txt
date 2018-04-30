cmake_minimum_required(VERSION 2.8)
project(GUM)
message("Building GUM in ${CMAKE_CURRENT_SOURCE_DIR}")

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/cmake)

# Check for Mathematica
include(cmake/FindMathematica.cmake)
if(Mathematica_FOUND)
  message("   Found Mathematica")
  if(Mathematica_WSTP_FOUND)
    message("   Found Wolfram Symbolic Transfer Protocol. Mathematica interface enabled.")
    set(HAVE_MATHEMATICA 1)
    set(MATHEMATICA_WSTP_H "${Mathematica_WSTP_INCLUDE_DIR}/wstp.h")
    set(MATHEMATICA_KERNEL "${Mathematica_KERNEL_EXECUTABLE}")
  else()
    message("  WSTP not found. Please make sure it is installed before attempting to use Mathematica interface.")
    set(HAVE_MATHEMATICA 0)
  endif()
else()
  message("   Mathematica not found. GUM won't work.")
  set(HAVE_MATHEMATICA 0)
endif()

if(Mathematica_FOUND AND Mathematica_WSTP_FOUND)
  set(LIBRARIES ${LIBRARIES} ${Mathematica_WSTP_LIBRARIES} ${Mathematica_MathLink_LIBRARY})
endif()

# Get FeynRules path
set(FR "" CACHE STRING "The path to the FeynRules module.")
get_filename_component(FEYNRULES_PATH "${FR}" PATH)

if(FEYNRULES_PATH)
  if(EXISTS "${FR}/FeynRules.m")
    message("   Found FeynRules.")
    set(FEYNRULES_PATH "${FR}")
    set(HAVE_FEYNRULES 1)
  else()
    set(HAVE_FEYNRULES 0)
    message(FATAL_ERROR "\n  FeynRules.m not found where you said it would be...")
  endif()
endif()  

# Get SARAH path
set(SARAH "" CACHE STRING "The path to the SARAH module.")
get_filename_component(SARAH_PATH "${SARAH}" PATH)

if(SARAH_PATH)
  if(EXISTS "${SARAH}/SARAH.m")
    message("   Found SARAH.")
    set(SARAH_PATH "${SARAH}")
    set(HAVE_SARAH 1)
  else()
    set(HAVE_SARAH 0)
    message(FATAL_ERROR "\n SARAH.m not found where you said it would be...")
  endif()
endif()

if(NOT HAVE_FEYNRULES AND NOT HAVE_SARAH)
  message(FATAL_ERROR "\n  You need to link GUM to either FeynRules or SARAH to get it to do stuff! \n  Use -D FR=path_to_feynrules or -D SARAH=path_to_sarah. (Or both.)")
endif()
  
# Configure cmake variables file
set(outhdr "${PROJECT_SOURCE_DIR}/include/cmake_variables.hpp")
configure_file("${PROJECT_SOURCE_DIR}/include/cmake_variables.hpp.in" ${outhdr})

include_directories(src include ${Mathematica_WSTP_INCLUDE_DIR})

# set source and header files used in all executables
set(HEADER_FILES  include/cmake_variables.hpp include/options.hpp include/feynrules.hpp include/sarah.hpp)

set(SRC_FILES src/mathematica.cpp)

set(SRC_GUM src/mathematica.cpp src/options.cpp src/feynrules.cpp src/sarah.cpp)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
set(CMAKE_CXX_FLAGS -L.)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")

add_executable(mathematica ${SRC_GUM} ${SRC_FILES_GEN} ${HEADER_FILES_GEN} ${HEADER_MB} ${TSIL_HEADER_FILE} ${SRC_FILES} ${HEADER_FILES} ${GENERATED_SRC_FILES})

TARGET_LINK_LIBRARIES(mathematica ${LIBRARIES})

if (CMAKE_COMPILER_IS_GNUCXX)
    message("Setting GCC flags")

    # Remove the -g flag when not debugging
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O2 -Wall -Wextra -W -pedantic -std=c++0x -O3")
else()
    message(status ": Setting MSVC flags")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /EHc-")
endif()
install(TARGETS mathematica DESTINATION ${CMAKE_BINARY_DIR}/lib)
