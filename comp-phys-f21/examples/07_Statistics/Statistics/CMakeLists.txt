# CMake file for CompPhys executable
cmake_minimum_required(VERSION 2.8)
project(projChi2Prob)

# Executable Name and Source Files
set( MY_EXE chi2prob )
set( MY_SRCS ${MY_EXE}.cpp )

# Libraries
set( CMAKE_LIBRARY_DIR $ENV{CPCODE}/library/lib )
set( CMAKE_LIBRARIES cpFunctions.so cpStatistics.so )
link_directories( ${CMAKE_LIBRARY_DIR} )

set( CMAKE_INCLUDE_DIR $ENV{CPCODE}/library/include )
include_directories( ${CMAKE_INCLUDE_DIR} )

# Generate executable and link to library
add_executable( ${MY_EXE} ${MY_SRCS} )
target_link_libraries( ${MY_EXE} ${CMAKE_LIBRARIES} )
