# CMAKE generated file: DO NOT EDIT!
# Generated by "MinGW Makefiles" Generator, CMake Version 3.12

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

SHELL = cmd.exe

# The CMake executable.
CMAKE_COMMAND = "C:\Program Files\JetBrains\CLion 2018.2.5\bin\cmake\win\bin\cmake.exe"

# The command to remove a file.
RM = "C:\Program Files\JetBrains\CLion 2018.2.5\bin\cmake\win\bin\cmake.exe" -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = C:\Users\Administrator\CLionProjects\Kwanty

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = C:\Users\Administrator\CLionProjects\Kwanty\cmake-build-debug

# Include any dependencies generated for this target.
include CMakeFiles/Kwanty.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/Kwanty.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/Kwanty.dir/flags.make

CMakeFiles/Kwanty.dir/main.cpp.obj: CMakeFiles/Kwanty.dir/flags.make
CMakeFiles/Kwanty.dir/main.cpp.obj: ../main.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=C:\Users\Administrator\CLionProjects\Kwanty\cmake-build-debug\CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/Kwanty.dir/main.cpp.obj"
	C:\MinGW\bin\g++.exe  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles\Kwanty.dir\main.cpp.obj -c C:\Users\Administrator\CLionProjects\Kwanty\main.cpp

CMakeFiles/Kwanty.dir/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/Kwanty.dir/main.cpp.i"
	C:\MinGW\bin\g++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E C:\Users\Administrator\CLionProjects\Kwanty\main.cpp > CMakeFiles\Kwanty.dir\main.cpp.i

CMakeFiles/Kwanty.dir/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/Kwanty.dir/main.cpp.s"
	C:\MinGW\bin\g++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S C:\Users\Administrator\CLionProjects\Kwanty\main.cpp -o CMakeFiles\Kwanty.dir\main.cpp.s

# Object files for target Kwanty
Kwanty_OBJECTS = \
"CMakeFiles/Kwanty.dir/main.cpp.obj"

# External object files for target Kwanty
Kwanty_EXTERNAL_OBJECTS =

Kwanty.exe: CMakeFiles/Kwanty.dir/main.cpp.obj
Kwanty.exe: CMakeFiles/Kwanty.dir/build.make
Kwanty.exe: CMakeFiles/Kwanty.dir/linklibs.rsp
Kwanty.exe: CMakeFiles/Kwanty.dir/objects1.rsp
Kwanty.exe: CMakeFiles/Kwanty.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=C:\Users\Administrator\CLionProjects\Kwanty\cmake-build-debug\CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable Kwanty.exe"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles\Kwanty.dir\link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/Kwanty.dir/build: Kwanty.exe

.PHONY : CMakeFiles/Kwanty.dir/build

CMakeFiles/Kwanty.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles\Kwanty.dir\cmake_clean.cmake
.PHONY : CMakeFiles/Kwanty.dir/clean

CMakeFiles/Kwanty.dir/depend:
	$(CMAKE_COMMAND) -E cmake_depends "MinGW Makefiles" C:\Users\Administrator\CLionProjects\Kwanty C:\Users\Administrator\CLionProjects\Kwanty C:\Users\Administrator\CLionProjects\Kwanty\cmake-build-debug C:\Users\Administrator\CLionProjects\Kwanty\cmake-build-debug C:\Users\Administrator\CLionProjects\Kwanty\cmake-build-debug\CMakeFiles\Kwanty.dir\DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/Kwanty.dir/depend

