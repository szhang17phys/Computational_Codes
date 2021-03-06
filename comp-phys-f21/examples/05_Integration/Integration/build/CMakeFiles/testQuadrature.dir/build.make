# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.18

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Disable VCS-based implicit rules.
% : %,v


# Disable VCS-based implicit rules.
% : RCS/%


# Disable VCS-based implicit rules.
% : RCS/%,v


# Disable VCS-based implicit rules.
% : SCCS/s.%


# Disable VCS-based implicit rules.
% : s.%


.SUFFIXES: .hpux_make_needs_suffix_list


# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /geode2/soft/hps/rhel7/cmake/gnu/3.18.4/bin/cmake

# The command to remove a file.
RM = /geode2/soft/hps/rhel7/cmake/gnu/3.18.4/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /N/u/szh2/Carbonate/CompPhys_szhang/comp-phys-f21/examples/05_Integration/Integration

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /N/u/szh2/Carbonate/CompPhys_szhang/comp-phys-f21/examples/05_Integration/Integration/build

# Include any dependencies generated for this target.
include CMakeFiles/testQuadrature.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/testQuadrature.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/testQuadrature.dir/flags.make

CMakeFiles/testQuadrature.dir/testQuadrature.cpp.o: CMakeFiles/testQuadrature.dir/flags.make
CMakeFiles/testQuadrature.dir/testQuadrature.cpp.o: ../testQuadrature.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/N/u/szh2/Carbonate/CompPhys_szhang/comp-phys-f21/examples/05_Integration/Integration/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/testQuadrature.dir/testQuadrature.cpp.o"
	/N/soft/rhel7/gcc/6.3.0/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/testQuadrature.dir/testQuadrature.cpp.o -c /N/u/szh2/Carbonate/CompPhys_szhang/comp-phys-f21/examples/05_Integration/Integration/testQuadrature.cpp

CMakeFiles/testQuadrature.dir/testQuadrature.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/testQuadrature.dir/testQuadrature.cpp.i"
	/N/soft/rhel7/gcc/6.3.0/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /N/u/szh2/Carbonate/CompPhys_szhang/comp-phys-f21/examples/05_Integration/Integration/testQuadrature.cpp > CMakeFiles/testQuadrature.dir/testQuadrature.cpp.i

CMakeFiles/testQuadrature.dir/testQuadrature.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/testQuadrature.dir/testQuadrature.cpp.s"
	/N/soft/rhel7/gcc/6.3.0/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /N/u/szh2/Carbonate/CompPhys_szhang/comp-phys-f21/examples/05_Integration/Integration/testQuadrature.cpp -o CMakeFiles/testQuadrature.dir/testQuadrature.cpp.s

# Object files for target testQuadrature
testQuadrature_OBJECTS = \
"CMakeFiles/testQuadrature.dir/testQuadrature.cpp.o"

# External object files for target testQuadrature
testQuadrature_EXTERNAL_OBJECTS =

testQuadrature: CMakeFiles/testQuadrature.dir/testQuadrature.cpp.o
testQuadrature: CMakeFiles/testQuadrature.dir/build.make
testQuadrature: CMakeFiles/testQuadrature.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/N/u/szh2/Carbonate/CompPhys_szhang/comp-phys-f21/examples/05_Integration/Integration/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable testQuadrature"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/testQuadrature.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/testQuadrature.dir/build: testQuadrature

.PHONY : CMakeFiles/testQuadrature.dir/build

CMakeFiles/testQuadrature.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/testQuadrature.dir/cmake_clean.cmake
.PHONY : CMakeFiles/testQuadrature.dir/clean

CMakeFiles/testQuadrature.dir/depend:
	cd /N/u/szh2/Carbonate/CompPhys_szhang/comp-phys-f21/examples/05_Integration/Integration/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /N/u/szh2/Carbonate/CompPhys_szhang/comp-phys-f21/examples/05_Integration/Integration /N/u/szh2/Carbonate/CompPhys_szhang/comp-phys-f21/examples/05_Integration/Integration /N/u/szh2/Carbonate/CompPhys_szhang/comp-phys-f21/examples/05_Integration/Integration/build /N/u/szh2/Carbonate/CompPhys_szhang/comp-phys-f21/examples/05_Integration/Integration/build /N/u/szh2/Carbonate/CompPhys_szhang/comp-phys-f21/examples/05_Integration/Integration/build/CMakeFiles/testQuadrature.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/testQuadrature.dir/depend

