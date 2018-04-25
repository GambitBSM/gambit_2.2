cmake_minimum_required(VERSION 2.8)
project(mass_builder)
message(status "${CMAKE_CURRENT_SOURCE_DIR}")

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/cmake)

# Check for Mathematica
include(cmake/FindMathematica.cmake)
if(Mathematica_FOUND)
  message("${BoldYellow}   Found Mathematica")
  if(Mathematica_WSTP_FOUND)
    message("${BoldYellow}   Found Wolfram Symbolic Transfer Protocol. Mathematica backends enabled")
    set(HAVE_MATHEMATICA 1)
    set(MATHEMATICA_WSTP_H "${Mathematica_WSTP_INCLUDE_DIR}/wstp.h")
    set(MATHEMATICA_KERNEL "${Mathematica_KERNEL_EXECUTABLE}")
  else()
    message("${BoldRed}  WSTP not found. Please make sure it is installed before attempting to use Mathematica backends")
    set(HAVE_MATHEMATICA 0)
  endif()
else()
  message("${BoldRed}   Mathematica not found. Backends using Mathematica will be disabled")
  set(HAVE_MATHEMATICA 0)
endif()

if(Mathematica_FOUND AND Mathematica_WSTP_FOUND)
  set(LIBRARIES ${LIBRARIES} ${Mathematica_WSTP_LIBRARIES} ${Mathematica_MathLink_LIBRARY})
endif()

# Configure cmake variables file
set(outhdr "${PROJECT_SOURCE_DIR}/include/cmake_variables.hpp")
configure_file("${PROJECT_SOURCE_DIR}/include/cmake_variables.hpp.in" ${outhdr})

include_directories(src include ${Mathematica_WSTP_INCLUDE_DIR})

# set source and header files used in all executables
set(HEADER_FILES  include/cmake_variables.hpp include/feynrules.hpp)

set(SRC_FILES src/mathematica.cpp)

set(SRC_GUM src/mathematica.cpp src/feynrules.cpp)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
set(CMAKE_CXX_FLAGS -L.)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")

add_executable(mathematica ${SRC_GUM} ${SRC_FILES_GEN} ${HEADER_FILES_GEN} ${HEADER_MB} ${TSIL_HEADER_FILE} ${SRC_FILES} ${HEADER_FILES} ${GENERATED_SRC_FILES})

TARGET_LINK_LIBRARIES(mathematica ${LIBRARIES})

if (CMAKE_COMPILER_IS_GNUCXX)
    message(status "Setting GCC flags")

    # Remove the -g flag when not debugging
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O2 -Wall -Wextra -W -pedantic -std=c++0x -O3")
else()
    message(status "Setting MSVC flags")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /EHc-")
endif()
install(TARGETS mathematica DESTINATION ${CMAKE_BINARY_DIR}/lib)
