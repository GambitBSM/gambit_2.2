cmake_minimum_required(VERSION 2.8)
project(GUM)
message("Building GUM in ${CMAKE_CURRENT_SOURCE_DIR}")

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/cmake)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/lib)

message("${CMAKE_CXX_COMPILER}")

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

# Python
FIND_PACKAGE(PythonInterp 2 REQUIRED)
if (PYTHON_VERSION_MAJOR GREATER 2)
  message(FATAL_ERROR "\n  GUM only supports Python 2.  \nIf you need to set the path to the Python interpreter manually, "
                      "please use -DPYTHON_EXECUTABLE:FILEPATH=path/to/preferred/python.")
endif()
FIND_PACKAGE(PythonLibs REQUIRED)
INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_DIRS})

# Add -fPIC for 64 bit systems
if(CMAKE_SYSTEM_PROCESSOR MATCHES "x86_64")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
endif()

# Check for C++11 and C++14 support
include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++14" COMPILER_SUPPORTS_CXX14)
if(COMPILER_SUPPORTS_CXX14)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")
else()
  CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
  if(COMPILER_SUPPORTS_CXX11)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
  else()
    CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)
    if(COMPILER_SUPPORTS_CXX0X)
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
    else()
      message(FATAL_ERROR "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler.")
    endif()
  endif()
endif()

# Check for Boost & dump FeynRules/SARAH shared object files. 
set(Boost_NO_BOOST_CMAKE ON)
FIND_PACKAGE(Boost 1.41)
if(Boost_FOUND)
  include_directories("${Boost_INCLUDE_DIR}")
  FIND_PACKAGE(Boost 1.41.0 COMPONENTS python)
  ADD_LIBRARY(fr SHARED src/feynrules.cpp)
  TARGET_LINK_LIBRARIES(fr ${LIBRARIES} ${PYTHON_LIBRARIES} ${Boost_LIBRARIES})
else()
  message(FATAL_ERROR "GUM requires Boost v1.41 or greater.\nPlease install a suitable version of Boost and rerun cmake.")
endif()
  
# Configure cmake variables file
set(outhdr "${PROJECT_SOURCE_DIR}/include/cmake_variables.hpp")
configure_file("${PROJECT_SOURCE_DIR}/include/cmake_variables.hpp.in" ${outhdr})

include_directories(src include ${Mathematica_WSTP_INCLUDE_DIR})

# set source and header files used in all executables
set(HEADER_FILES include/cmake_variables.hpp include/options.hpp include/feynrules.hpp include/sarah.hpp)
set(SRC_FILES src/mathematica.cpp)
set(SRC_GUM src/mathematica.cpp src/options.cpp src/feynrules.cpp src/sarah.cpp)


set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
set(CMAKE_CXX_FLAGS -L.)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")

add_executable(mathematica ${SRC_GUM} ${SRC_FILES_GEN} ${HEADER_FILES_GEN} ${HEADER_MB} ${SRC_FILES} ${HEADER_FILES} ${GENERATED_SRC_FILES})

TARGET_LINK_LIBRARIES(mathematica ${LIBRARIES} ${PYTHON_LIBRARIES} ${Boost_LIBRARIES})

if (CMAKE_COMPILER_IS_GNUCXX)
    message("Setting GCC flags")

    # Remove the -g flag when not debugging
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O2 -Wall -Wextra -W -pedantic -std=c++0x -O3")
else()
    message(status ": Setting MSVC flags")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /EHc-")
endif()
install(TARGETS mathematica DESTINATION ${CMAKE_BINARY_DIR}/lib)
