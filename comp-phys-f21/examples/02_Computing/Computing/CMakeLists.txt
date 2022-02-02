# CMake file for CompPhys executable
cmake_minimum_required(VERSION 2.8)
project(projComputing)

# Executable Name and Source Files
set( MY_EXE testComputing )
set( MY_SRCS ${MY_EXE}.cpp )

# Generate executable and link to library
add_executable( ${MY_EXE} ${MY_SRCS} )
