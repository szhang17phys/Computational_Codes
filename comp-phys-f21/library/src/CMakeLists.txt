# CMake file for CompPhys libraries
cmake_minimum_required(VERSION 2.8)
project(cplibraries)
set(CMAKE_BUILD_TYPE Release)

# Directory structure: store libraries to lib/ directory
set(CMAKE_BINARY_DIR ${CMAKE_SOURCE_DIR}/../lib) 
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})
set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR})
 
# Header files directory
include_directories(../include)

# List of Libraries to make
set(LIBRARIES cpDoubleRep cpODEmodel cpODEstep cpInterpolate cpFunctions cpQuadrature cpRandom cpStatistics cpMatrix cpFourier cpPDE cpGenSample )

# Generate the shared libraries from the sources
foreach(library ${LIBRARIES})
  add_library(${library} SHARED ${library}.cpp)
endforeach(library)
