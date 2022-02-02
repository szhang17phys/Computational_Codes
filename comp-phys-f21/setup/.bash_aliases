#!/bin/sh
#---------------------------------------------------------------------------------------------
# Aliases for CompPhys bash shell
# Modifications
# 27-Aug-2018  add emacs and root modules
# 05-Jun-2018  created
#---------------------------------------------------------------------------------------------

# Aliases
alias purge="rm *~; rm .*~"

# CMake/Eclipse aliases, to run:
#  $ cmake_eclipse <path to CMakeLists.txt directory>
alias cmake_eclipse="cmake -G\"Eclipse CDT4 - Unix Makefiles\" -D CMAKE_BUILD_TYPE=Debug -D CMAKE_ECLIPSE_VERSION=4.7"

# Environment Variables
export CPHOME=$HOME/CompPhys
export CPCODE=$HOME/CompPhys/comp-phys-f21

# Modules
module load emacs
module load root
