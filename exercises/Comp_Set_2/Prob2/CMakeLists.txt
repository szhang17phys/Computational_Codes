# CMake file for CompPhys executable
cmake_minimum_required(VERSION 2.8)
project(projQuadrature)

# Executable Name and Source Files
set( MY_EXE prob_2 )
set( MY_SRCS ${MY_EXE}.cpp )

#Find ROOT library; to use root! --szhang, Nov 1, 2021---
find_package(ROOT REQUIRED COMPONENTS Core Hist Gpad)

# Libraries
set( CMAKE_LIBRARY_DIR $ENV{CPCODE}/library/lib )
set( CMAKE_LIBRARIES cpQuadrature.so cpFunctions.so cpInterpolate.so )
link_directories( ${CMAKE_LIBRARY_DIR} )

set( CMAKE_INCLUDE_DIR $ENV{CPCODE}/library/include )
include_directories( ${CMAKE_INCLUDE_DIR} )

# Generate executable and link to library
add_executable( ${MY_EXE} ${MY_SRCS} )
target_link_libraries( ${MY_EXE} ${CMAKE_LIBRARIES} ROOT::Core ROOT::Hist ROOT::Gpad )# ROOT:... is about root -szhang, Nov 1, 2021---
