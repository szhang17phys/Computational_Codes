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
CMAKE_SOURCE_DIR = /N/u/szh2/Carbonate/CompPhys_szhang/comp-phys-f21/library/src

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /N/u/szh2/Carbonate/CompPhys_szhang/comp-phys-f21/library/build

# Include any dependencies generated for this target.
include CMakeFiles/cpGenSample.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/cpGenSample.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/cpGenSample.dir/flags.make

CMakeFiles/cpGenSample.dir/cpGenSample.cpp.o: CMakeFiles/cpGenSample.dir/flags.make
CMakeFiles/cpGenSample.dir/cpGenSample.cpp.o: /N/u/szh2/Carbonate/CompPhys_szhang/comp-phys-f21/library/src/cpGenSample.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/N/u/szh2/Carbonate/CompPhys_szhang/comp-phys-f21/library/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/cpGenSample.dir/cpGenSample.cpp.o"
	/N/soft/rhel7/gcc/6.3.0/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/cpGenSample.dir/cpGenSample.cpp.o -c /N/u/szh2/Carbonate/CompPhys_szhang/comp-phys-f21/library/src/cpGenSample.cpp

CMakeFiles/cpGenSample.dir/cpGenSample.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/cpGenSample.dir/cpGenSample.cpp.i"
	/N/soft/rhel7/gcc/6.3.0/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /N/u/szh2/Carbonate/CompPhys_szhang/comp-phys-f21/library/src/cpGenSample.cpp > CMakeFiles/cpGenSample.dir/cpGenSample.cpp.i

CMakeFiles/cpGenSample.dir/cpGenSample.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/cpGenSample.dir/cpGenSample.cpp.s"
	/N/soft/rhel7/gcc/6.3.0/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /N/u/szh2/Carbonate/CompPhys_szhang/comp-phys-f21/library/src/cpGenSample.cpp -o CMakeFiles/cpGenSample.dir/cpGenSample.cpp.s

# Object files for target cpGenSample
cpGenSample_OBJECTS = \
"CMakeFiles/cpGenSample.dir/cpGenSample.cpp.o"

# External object files for target cpGenSample
cpGenSample_EXTERNAL_OBJECTS =

/N/u/szh2/Carbonate/CompPhys_szhang/comp-phys-f21/library/lib/libcpGenSample.so: CMakeFiles/cpGenSample.dir/cpGenSample.cpp.o
/N/u/szh2/Carbonate/CompPhys_szhang/comp-phys-f21/library/lib/libcpGenSample.so: CMakeFiles/cpGenSample.dir/build.make
/N/u/szh2/Carbonate/CompPhys_szhang/comp-phys-f21/library/lib/libcpGenSample.so: CMakeFiles/cpGenSample.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/N/u/szh2/Carbonate/CompPhys_szhang/comp-phys-f21/library/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX shared library /N/u/szh2/Carbonate/CompPhys_szhang/comp-phys-f21/library/lib/libcpGenSample.so"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/cpGenSample.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/cpGenSample.dir/build: /N/u/szh2/Carbonate/CompPhys_szhang/comp-phys-f21/library/lib/libcpGenSample.so

.PHONY : CMakeFiles/cpGenSample.dir/build

CMakeFiles/cpGenSample.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/cpGenSample.dir/cmake_clean.cmake
.PHONY : CMakeFiles/cpGenSample.dir/clean

CMakeFiles/cpGenSample.dir/depend:
	cd /N/u/szh2/Carbonate/CompPhys_szhang/comp-phys-f21/library/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /N/u/szh2/Carbonate/CompPhys_szhang/comp-phys-f21/library/src /N/u/szh2/Carbonate/CompPhys_szhang/comp-phys-f21/library/src /N/u/szh2/Carbonate/CompPhys_szhang/comp-phys-f21/library/build /N/u/szh2/Carbonate/CompPhys_szhang/comp-phys-f21/library/build /N/u/szh2/Carbonate/CompPhys_szhang/comp-phys-f21/library/build/CMakeFiles/cpGenSample.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/cpGenSample.dir/depend

