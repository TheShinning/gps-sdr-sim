# CMAKE generated file: DO NOT EDIT!
# Generated by "MinGW Makefiles" Generator, CMake Version 3.23

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

SHELL = cmd.exe

# The CMake executable.
CMAKE_COMMAND = E:\software\codesoftware\Cmake\bin\cmake.exe

# The command to remove a file.
RM = E:\software\codesoftware\Cmake\bin\cmake.exe -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = H:\winmove\project\pppar\pppartest\GARPS

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = H:\winmove\project\pppar\pppartest\GARPS\build

# Include any dependencies generated for this target.
include libQc/CMakeFiles/libQc.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include libQc/CMakeFiles/libQc.dir/compiler_depend.make

# Include the progress variables for this target.
include libQc/CMakeFiles/libQc.dir/progress.make

# Include the compile flags for this target's objects.
include libQc/CMakeFiles/libQc.dir/flags.make

libQc/CMakeFiles/libQc.dir/qc.cpp.obj: libQc/CMakeFiles/libQc.dir/flags.make
libQc/CMakeFiles/libQc.dir/qc.cpp.obj: libQc/CMakeFiles/libQc.dir/includes_CXX.rsp
libQc/CMakeFiles/libQc.dir/qc.cpp.obj: ../src/LibQc/qc.cpp
libQc/CMakeFiles/libQc.dir/qc.cpp.obj: libQc/CMakeFiles/libQc.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=H:\winmove\project\pppar\pppartest\GARPS\build\CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object libQc/CMakeFiles/libQc.dir/qc.cpp.obj"
	cd /d H:\winmove\project\pppar\pppartest\GARPS\build\libQc && E:\software\codesoftware\mingw64\bin\g++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT libQc/CMakeFiles/libQc.dir/qc.cpp.obj -MF CMakeFiles\libQc.dir\qc.cpp.obj.d -o CMakeFiles\libQc.dir\qc.cpp.obj -c H:\winmove\project\pppar\pppartest\GARPS\src\LibQc\qc.cpp

libQc/CMakeFiles/libQc.dir/qc.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/libQc.dir/qc.cpp.i"
	cd /d H:\winmove\project\pppar\pppartest\GARPS\build\libQc && E:\software\codesoftware\mingw64\bin\g++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E H:\winmove\project\pppar\pppartest\GARPS\src\LibQc\qc.cpp > CMakeFiles\libQc.dir\qc.cpp.i

libQc/CMakeFiles/libQc.dir/qc.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/libQc.dir/qc.cpp.s"
	cd /d H:\winmove\project\pppar\pppartest\GARPS\build\libQc && E:\software\codesoftware\mingw64\bin\g++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S H:\winmove\project\pppar\pppartest\GARPS\src\LibQc\qc.cpp -o CMakeFiles\libQc.dir\qc.cpp.s

# Object files for target libQc
libQc_OBJECTS = \
"CMakeFiles/libQc.dir/qc.cpp.obj"

# External object files for target libQc
libQc_EXTERNAL_OBJECTS =

Lib/liblibQcd.a: libQc/CMakeFiles/libQc.dir/qc.cpp.obj
Lib/liblibQcd.a: libQc/CMakeFiles/libQc.dir/build.make
Lib/liblibQcd.a: libQc/CMakeFiles/libQc.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=H:\winmove\project\pppar\pppartest\GARPS\build\CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX static library ..\Lib\liblibQcd.a"
	cd /d H:\winmove\project\pppar\pppartest\GARPS\build\libQc && $(CMAKE_COMMAND) -P CMakeFiles\libQc.dir\cmake_clean_target.cmake
	cd /d H:\winmove\project\pppar\pppartest\GARPS\build\libQc && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles\libQc.dir\link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
libQc/CMakeFiles/libQc.dir/build: Lib/liblibQcd.a
.PHONY : libQc/CMakeFiles/libQc.dir/build

libQc/CMakeFiles/libQc.dir/clean:
	cd /d H:\winmove\project\pppar\pppartest\GARPS\build\libQc && $(CMAKE_COMMAND) -P CMakeFiles\libQc.dir\cmake_clean.cmake
.PHONY : libQc/CMakeFiles/libQc.dir/clean

libQc/CMakeFiles/libQc.dir/depend:
	$(CMAKE_COMMAND) -E cmake_depends "MinGW Makefiles" H:\winmove\project\pppar\pppartest\GARPS H:\winmove\project\pppar\pppartest\GARPS\src\LibQc H:\winmove\project\pppar\pppartest\GARPS\build H:\winmove\project\pppar\pppartest\GARPS\build\libQc H:\winmove\project\pppar\pppartest\GARPS\build\libQc\CMakeFiles\libQc.dir\DependInfo.cmake --color=$(COLOR)
.PHONY : libQc/CMakeFiles/libQc.dir/depend

