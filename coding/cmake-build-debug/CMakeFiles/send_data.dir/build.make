# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.10

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = "/mnt/c/Users/1513 MXTI/ITA/PROF/2SEM/CTC-34-L/coding"

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = "/mnt/c/Users/1513 MXTI/ITA/PROF/2SEM/CTC-34-L/coding/cmake-build-debug"

# Include any dependencies generated for this target.
include CMakeFiles/send_data.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/send_data.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/send_data.dir/flags.make

CMakeFiles/send_data.dir/send_data.cpp.o: CMakeFiles/send_data.dir/flags.make
CMakeFiles/send_data.dir/send_data.cpp.o: ../send_data.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir="/mnt/c/Users/1513 MXTI/ITA/PROF/2SEM/CTC-34-L/coding/cmake-build-debug/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/send_data.dir/send_data.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/send_data.dir/send_data.cpp.o -c "/mnt/c/Users/1513 MXTI/ITA/PROF/2SEM/CTC-34-L/coding/send_data.cpp"

CMakeFiles/send_data.dir/send_data.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/send_data.dir/send_data.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E "/mnt/c/Users/1513 MXTI/ITA/PROF/2SEM/CTC-34-L/coding/send_data.cpp" > CMakeFiles/send_data.dir/send_data.cpp.i

CMakeFiles/send_data.dir/send_data.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/send_data.dir/send_data.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S "/mnt/c/Users/1513 MXTI/ITA/PROF/2SEM/CTC-34-L/coding/send_data.cpp" -o CMakeFiles/send_data.dir/send_data.cpp.s

CMakeFiles/send_data.dir/send_data.cpp.o.requires:

.PHONY : CMakeFiles/send_data.dir/send_data.cpp.o.requires

CMakeFiles/send_data.dir/send_data.cpp.o.provides: CMakeFiles/send_data.dir/send_data.cpp.o.requires
	$(MAKE) -f CMakeFiles/send_data.dir/build.make CMakeFiles/send_data.dir/send_data.cpp.o.provides.build
.PHONY : CMakeFiles/send_data.dir/send_data.cpp.o.provides

CMakeFiles/send_data.dir/send_data.cpp.o.provides.build: CMakeFiles/send_data.dir/send_data.cpp.o


# Object files for target send_data
send_data_OBJECTS = \
"CMakeFiles/send_data.dir/send_data.cpp.o"

# External object files for target send_data
send_data_EXTERNAL_OBJECTS =

send_data: CMakeFiles/send_data.dir/send_data.cpp.o
send_data: CMakeFiles/send_data.dir/build.make
send_data: CMakeFiles/send_data.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir="/mnt/c/Users/1513 MXTI/ITA/PROF/2SEM/CTC-34-L/coding/cmake-build-debug/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable send_data"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/send_data.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/send_data.dir/build: send_data

.PHONY : CMakeFiles/send_data.dir/build

CMakeFiles/send_data.dir/requires: CMakeFiles/send_data.dir/send_data.cpp.o.requires

.PHONY : CMakeFiles/send_data.dir/requires

CMakeFiles/send_data.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/send_data.dir/cmake_clean.cmake
.PHONY : CMakeFiles/send_data.dir/clean

CMakeFiles/send_data.dir/depend:
	cd "/mnt/c/Users/1513 MXTI/ITA/PROF/2SEM/CTC-34-L/coding/cmake-build-debug" && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" "/mnt/c/Users/1513 MXTI/ITA/PROF/2SEM/CTC-34-L/coding" "/mnt/c/Users/1513 MXTI/ITA/PROF/2SEM/CTC-34-L/coding" "/mnt/c/Users/1513 MXTI/ITA/PROF/2SEM/CTC-34-L/coding/cmake-build-debug" "/mnt/c/Users/1513 MXTI/ITA/PROF/2SEM/CTC-34-L/coding/cmake-build-debug" "/mnt/c/Users/1513 MXTI/ITA/PROF/2SEM/CTC-34-L/coding/cmake-build-debug/CMakeFiles/send_data.dir/DependInfo.cmake" --color=$(COLOR)
.PHONY : CMakeFiles/send_data.dir/depend
