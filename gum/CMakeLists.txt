cmake_minimum_required(VERSION 2.8)
project(GUM)
message("Building GUM in ${CMAKE_CURRENT_SOURCE_DIR}")
include(ExternalProject)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/cmake)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/lib)

# Define the main GUM directory
set(GUM_DIR ${PROJECT_SOURCE_DIR})

# Define the sed command to use differently for OSX and linux.  Delete this when properly integrated into GAMBIT cmake system, as it exists there already.
if (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
  set(dashi "-i ''")
else()
  set(dashi "-i")
endif()

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

# Download FeynRules
EXTERNALPROJECT_ADD(
    FeynRules
    URL http://feynrules.irmp.ucl.ac.be/downloads/feynrules-current.tar.gz
    UPDATE_COMMAND ""
    PATCH_COMMAND ""
    SOURCE_DIR ${CMAKE_SOURCE_DIR}/contrib/FeynRules
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND ""
#    INSTALL_COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_SOURCE_DIR}/Models ${CMAKE_SOURCE_DIR}/contrib/FeynRules/Models
)

# Download SARAH
set(SARAH_VERSION 4.14.0)
EXTERNALPROJECT_ADD(
    SARAH
    URL https://sarah.hepforge.org/downloads/?f=SARAH-${SARAH_VERSION}.tar.gz
    URL_MD5 850b74625e531b93fd43a32c181c825b
    UPDATE_COMMAND ""
    PATCH_COMMAND ""
    SOURCE_DIR ${CMAKE_SOURCE_DIR}/contrib/SARAH
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND ""
)

# Download MadGraph
EXTERNALPROJECT_ADD(
    MadGraph
    URL https://launchpad.net/mg5amcnlo/2.0/2.6.x/+download/MG5_aMC_v2.6.7.tar.gz
    UPDATE_COMMAND ""
    PATCH_COMMAND sed ${dashi} -e "s/# eps_viewer = None/eps_viewer = true/g" input/mg5_configuration.txt
    SOURCE_DIR ${CMAKE_SOURCE_DIR}/contrib/MadGraph
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND ""
)

# Download Pythia
set(BASE_PYTHIA_VERSION 212)
set(PYTHIA_MD5 0886d1b2827d8f0cd2ae69b925045f40)
EXTERNALPROJECT_ADD(
    Pythia
    URL http://home.thep.lu.se/~torbjorn/pythia8/pythia8${BASE_PYTHIA_VERSION}.tgz
    URL_MD5 ${PYTHIA_MD5}
    UPDATE_COMMAND ""
    PATCH_COMMAND  ""
    SOURCE_DIR ${CMAKE_SOURCE_DIR}/contrib/pythia
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND ""
)

# Download SPheno
set(SPHENO_VERSION 4.0.3)
set(dl http://www.hepforge.org/archive/spheno/SPheno-${SPHENO_VERSION}.tar.gz)
set(md5 "64787d6c8ce03cac38aec53d34ac46ad")
set(SPHENO_DIR "${CMAKE_SOURCE_DIR}/contrib/SPheno")
ExternalProject_Add(
  SPheno
  URL ${dl}
  URL_MD5 ${md5}
  SOURCE_DIR ${SPHENO_DIR}
  PATCH_COMMAND ""
  CONFIGURE_COMMAND ""
  BUILD_COMMAND ""
  INSTALL_COMMAND ""
)

# Set paths
set(FEYNRULES_PATH "${CMAKE_SOURCE_DIR}/contrib/FeynRules/")
set(HAVE_FEYNRULES 1)
set(SARAH_PATH "${CMAKE_SOURCE_DIR}/contrib/SARAH/")
set(HAVE_SARAH 1)

# Python
FIND_PACKAGE(PythonInterp 2 REQUIRED)
if (PYTHON_VERSION_MAJOR GREATER 2)
  message(FATAL_ERROR "\n  GUM only supports Python 2.  \nIf you need to set the path to the Python interpreter manually, "
                      "please use -DPYTHON_EXECUTABLE:FILEPATH=path/to/preferred/python.")
endif()
FIND_PACKAGE(PythonLibs REQUIRED)
INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_DIRS})
message("Using Python interpreter version ${PYTHON_VERSION_STRING} for build.")

# Simple function to find specific Python modules
macro(find_python_module module)
  execute_process(COMMAND ${PYTHON_EXECUTABLE} -c "import ${module}" RESULT_VARIABLE return_value ERROR_QUIET)
  if (NOT return_value)
    message(STATUS "Found Python module ${module}.")
    set(PY_${module}_FOUND TRUE)
  else()
    if(${ARGC} GREATER 1 AND ${ARGV1} STREQUAL "REQUIRED")
      message(FATAL_ERROR "-- FAILED to find Python module ${module}.")
    else()
      message(STATUS "FAILED to find Python module ${module}.")
    endif()
  endif()
endmacro()

# Check for required Python libraries
foreach(module yaml os re datetime sys getopt shutil itertools collections copy numpy distutils filecmp glob six argparse)
  find_python_module(${module} REQUIRED)
endforeach()

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

# Check for Boost & drop FeynRules/SARAH shared object files.
# Note that for Boost 1.67 or later, we'll need to match the Python version with the component, as e.g. COMPONENTS "python36".
set(Boost_NO_BOOST_CMAKE ON)
find_package(Boost 1.41.0 COMPONENTS python filesystem)
if(Boost_FOUND)
  set( CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -Wl,--no-undefined" )
  include_directories("${Boost_INCLUDE_DIR}")
  add_library(fr SHARED src/feynrules.cpp src/math_package.cpp)
  target_link_libraries(fr ${LIBRARIES} ${PYTHON_LIBRARIES} ${Boost_LIBRARIES})
  add_library(sarah SHARED src/sarah.cpp src/math_package.cpp)
  target_link_libraries(sarah ${LIBRARIES} ${PYTHON_LIBRARIES} ${Boost_LIBRARIES})
else()
  message(FATAL_ERROR "GUM requires Boost v1.41 or greater.\nPlease install a suitable version of Boost and re-run cmake.")
endif()

# Configure cmake variables files
set(outfile "${PROJECT_SOURCE_DIR}/include/cmake_variables.hpp")
configure_file("${outfile}.in" ${outfile})
set(outfile "${PROJECT_SOURCE_DIR}/src/cmake_variables.py")
configure_file("${outfile}.in" ${outfile})

include_directories(src include ${Mathematica_WSTP_INCLUDE_DIR})

# set source and header files used in all executables
set(HEADER_FILES include/cmake_variables.hpp include/options.hpp include/math_package.hpp include/feynrules.hpp include/sarah.hpp)
set(SRC_FILES src/mathematica.cpp)
set(SRC_GUM src/mathematica.cpp src/options.cpp)
set(GUM_LIBS sarah fr)


set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
set(CMAKE_CXX_FLAGS -L.)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")

add_executable(mathematica ${SRC_GUM} ${HEADER_FILES})

TARGET_LINK_LIBRARIES(mathematica ${LIBRARIES} ${PYTHON_LIBRARIES} ${Boost_LIBRARIES} ${GUM_LIBS})

if (CMAKE_COMPILER_IS_GNUCXX)
    message("Setting GCC flags")

    # Remove the -g flag when not debugging
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O2 -Wall -Wextra -W -pedantic -std=c++0x -O3")
else()
    message(status ": Setting MSVC flags")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /EHc-")
endif()
install(TARGETS mathematica DESTINATION ${CMAKE_BINARY_DIR}/lib)
