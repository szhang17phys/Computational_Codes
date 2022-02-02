# CMake file for CompPhys executable
cmake_minimum_required(VERSION 2.8)
project(rooplot)

# Find ROOT libraries
find_package(ROOT REQUIRED COMPONENTS Core Hist Gpad)

# Executable Name and Source Files
set( MY_EXE plot )
set( MY_SRCS ${MY_EXE}.cpp )

# Generate executable and link to library
add_executable( ${MY_EXE} ${MY_SRCS} )
target_link_libraries( ${MY_EXE} ROOT::Core ROOT::Hist ROOT::Gpad )
