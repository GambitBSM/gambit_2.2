cmake_minimum_required(VERSION 2.8)
project(GUM)
message("Building GUM in ${CMAKE_CURRENT_SOURCE_DIR}")
include(ExternalProject)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/cmake)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/lib)

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
    URL_MD5 95c70262ac434d48346ef2f97165e76b
    
    UPDATE_COMMAND ""
    PATCH_COMMAND ""
    
    SOURCE_DIR "${CMAKE_SOURCE_DIR}/contrib/FeynRules"
    CONFIGURE_COMMAND ""
    
    BUILD_COMMAND ""
    
    INSTALL_COMMAND ""
    
    TEST_COMMAND ""
)

# Download SARAH
EXTERNALPROJECT_ADD(
    SARAH
    
    URL http://www.hepforge.org/archive/sarah/SARAH-4.12.3.tar.gz
    URL_MD5 be173d937e88f18df81dd29a08efc513
    
    UPDATE_COMMAND ""
    PATCH_COMMAND ""
    
    SOURCE_DIR "${CMAKE_SOURCE_DIR}/contrib/SARAH"
    CONFIGURE_COMMAND ""
    
    BUILD_COMMAND ""
    
    INSTALL_COMMAND ""
    
    TEST_COMMAND ""
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
  ADD_LIBRARY(sarah SHARED src/sarah.cpp)
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
