# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.16

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
CMAKE_COMMAND = /home/manuel/.local/share/JetBrains/Toolbox/apps/CLion/ch-0/201.7846.88/bin/cmake/linux/bin/cmake

# The command to remove a file.
RM = /home/manuel/.local/share/JetBrains/Toolbox/apps/CLion/ch-0/201.7846.88/bin/cmake/linux/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /media/sf_CS/git/koi_0x4249/SS20/WR/sc_heat_equation

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /media/sf_CS/git/koi_0x4249/SS20/WR/sc_heat_equation/cmake-build-debug

# Include any dependencies generated for this target.
include external/pmp/CMakeFiles/pmp.dir/depend.make

# Include the progress variables for this target.
include external/pmp/CMakeFiles/pmp.dir/progress.make

# Include the compile flags for this target's objects.
include external/pmp/CMakeFiles/pmp.dir/flags.make

external/pmp/CMakeFiles/pmp.dir/Shader.cpp.o: external/pmp/CMakeFiles/pmp.dir/flags.make
external/pmp/CMakeFiles/pmp.dir/Shader.cpp.o: ../external/pmp/Shader.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/media/sf_CS/git/koi_0x4249/SS20/WR/sc_heat_equation/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object external/pmp/CMakeFiles/pmp.dir/Shader.cpp.o"
	cd /media/sf_CS/git/koi_0x4249/SS20/WR/sc_heat_equation/cmake-build-debug/external/pmp && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/pmp.dir/Shader.cpp.o -c /media/sf_CS/git/koi_0x4249/SS20/WR/sc_heat_equation/external/pmp/Shader.cpp

external/pmp/CMakeFiles/pmp.dir/Shader.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/pmp.dir/Shader.cpp.i"
	cd /media/sf_CS/git/koi_0x4249/SS20/WR/sc_heat_equation/cmake-build-debug/external/pmp && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /media/sf_CS/git/koi_0x4249/SS20/WR/sc_heat_equation/external/pmp/Shader.cpp > CMakeFiles/pmp.dir/Shader.cpp.i

external/pmp/CMakeFiles/pmp.dir/Shader.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/pmp.dir/Shader.cpp.s"
	cd /media/sf_CS/git/koi_0x4249/SS20/WR/sc_heat_equation/cmake-build-debug/external/pmp && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /media/sf_CS/git/koi_0x4249/SS20/WR/sc_heat_equation/external/pmp/Shader.cpp -o CMakeFiles/pmp.dir/Shader.cpp.s

external/pmp/CMakeFiles/pmp.dir/Window.cpp.o: external/pmp/CMakeFiles/pmp.dir/flags.make
external/pmp/CMakeFiles/pmp.dir/Window.cpp.o: ../external/pmp/Window.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/media/sf_CS/git/koi_0x4249/SS20/WR/sc_heat_equation/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object external/pmp/CMakeFiles/pmp.dir/Window.cpp.o"
	cd /media/sf_CS/git/koi_0x4249/SS20/WR/sc_heat_equation/cmake-build-debug/external/pmp && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/pmp.dir/Window.cpp.o -c /media/sf_CS/git/koi_0x4249/SS20/WR/sc_heat_equation/external/pmp/Window.cpp

external/pmp/CMakeFiles/pmp.dir/Window.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/pmp.dir/Window.cpp.i"
	cd /media/sf_CS/git/koi_0x4249/SS20/WR/sc_heat_equation/cmake-build-debug/external/pmp && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /media/sf_CS/git/koi_0x4249/SS20/WR/sc_heat_equation/external/pmp/Window.cpp > CMakeFiles/pmp.dir/Window.cpp.i

external/pmp/CMakeFiles/pmp.dir/Window.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/pmp.dir/Window.cpp.s"
	cd /media/sf_CS/git/koi_0x4249/SS20/WR/sc_heat_equation/cmake-build-debug/external/pmp && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /media/sf_CS/git/koi_0x4249/SS20/WR/sc_heat_equation/external/pmp/Window.cpp -o CMakeFiles/pmp.dir/Window.cpp.s

# Object files for target pmp
pmp_OBJECTS = \
"CMakeFiles/pmp.dir/Shader.cpp.o" \
"CMakeFiles/pmp.dir/Window.cpp.o"

# External object files for target pmp
pmp_EXTERNAL_OBJECTS =

libpmp.a: external/pmp/CMakeFiles/pmp.dir/Shader.cpp.o
libpmp.a: external/pmp/CMakeFiles/pmp.dir/Window.cpp.o
libpmp.a: external/pmp/CMakeFiles/pmp.dir/build.make
libpmp.a: external/pmp/CMakeFiles/pmp.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/media/sf_CS/git/koi_0x4249/SS20/WR/sc_heat_equation/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Linking CXX static library ../../libpmp.a"
	cd /media/sf_CS/git/koi_0x4249/SS20/WR/sc_heat_equation/cmake-build-debug/external/pmp && $(CMAKE_COMMAND) -P CMakeFiles/pmp.dir/cmake_clean_target.cmake
	cd /media/sf_CS/git/koi_0x4249/SS20/WR/sc_heat_equation/cmake-build-debug/external/pmp && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/pmp.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
external/pmp/CMakeFiles/pmp.dir/build: libpmp.a

.PHONY : external/pmp/CMakeFiles/pmp.dir/build

external/pmp/CMakeFiles/pmp.dir/clean:
	cd /media/sf_CS/git/koi_0x4249/SS20/WR/sc_heat_equation/cmake-build-debug/external/pmp && $(CMAKE_COMMAND) -P CMakeFiles/pmp.dir/cmake_clean.cmake
.PHONY : external/pmp/CMakeFiles/pmp.dir/clean

external/pmp/CMakeFiles/pmp.dir/depend:
	cd /media/sf_CS/git/koi_0x4249/SS20/WR/sc_heat_equation/cmake-build-debug && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /media/sf_CS/git/koi_0x4249/SS20/WR/sc_heat_equation /media/sf_CS/git/koi_0x4249/SS20/WR/sc_heat_equation/external/pmp /media/sf_CS/git/koi_0x4249/SS20/WR/sc_heat_equation/cmake-build-debug /media/sf_CS/git/koi_0x4249/SS20/WR/sc_heat_equation/cmake-build-debug/external/pmp /media/sf_CS/git/koi_0x4249/SS20/WR/sc_heat_equation/cmake-build-debug/external/pmp/CMakeFiles/pmp.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : external/pmp/CMakeFiles/pmp.dir/depend
