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
include src/CMakeFiles/diffusion.dir/depend.make

# Include the progress variables for this target.
include src/CMakeFiles/diffusion.dir/progress.make

# Include the compile flags for this target's objects.
include src/CMakeFiles/diffusion.dir/flags.make

src/CMakeFiles/diffusion.dir/HeatEquationViewer.cpp.o: src/CMakeFiles/diffusion.dir/flags.make
src/CMakeFiles/diffusion.dir/HeatEquationViewer.cpp.o: ../src/HeatEquationViewer.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/media/sf_CS/git/koi_0x4249/SS20/WR/sc_heat_equation/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object src/CMakeFiles/diffusion.dir/HeatEquationViewer.cpp.o"
	cd /media/sf_CS/git/koi_0x4249/SS20/WR/sc_heat_equation/cmake-build-debug/src && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/diffusion.dir/HeatEquationViewer.cpp.o -c /media/sf_CS/git/koi_0x4249/SS20/WR/sc_heat_equation/src/HeatEquationViewer.cpp

src/CMakeFiles/diffusion.dir/HeatEquationViewer.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/diffusion.dir/HeatEquationViewer.cpp.i"
	cd /media/sf_CS/git/koi_0x4249/SS20/WR/sc_heat_equation/cmake-build-debug/src && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /media/sf_CS/git/koi_0x4249/SS20/WR/sc_heat_equation/src/HeatEquationViewer.cpp > CMakeFiles/diffusion.dir/HeatEquationViewer.cpp.i

src/CMakeFiles/diffusion.dir/HeatEquationViewer.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/diffusion.dir/HeatEquationViewer.cpp.s"
	cd /media/sf_CS/git/koi_0x4249/SS20/WR/sc_heat_equation/cmake-build-debug/src && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /media/sf_CS/git/koi_0x4249/SS20/WR/sc_heat_equation/src/HeatEquationViewer.cpp -o CMakeFiles/diffusion.dir/HeatEquationViewer.cpp.s

src/CMakeFiles/diffusion.dir/heat.cpp.o: src/CMakeFiles/diffusion.dir/flags.make
src/CMakeFiles/diffusion.dir/heat.cpp.o: ../src/heat.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/media/sf_CS/git/koi_0x4249/SS20/WR/sc_heat_equation/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object src/CMakeFiles/diffusion.dir/heat.cpp.o"
	cd /media/sf_CS/git/koi_0x4249/SS20/WR/sc_heat_equation/cmake-build-debug/src && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/diffusion.dir/heat.cpp.o -c /media/sf_CS/git/koi_0x4249/SS20/WR/sc_heat_equation/src/heat.cpp

src/CMakeFiles/diffusion.dir/heat.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/diffusion.dir/heat.cpp.i"
	cd /media/sf_CS/git/koi_0x4249/SS20/WR/sc_heat_equation/cmake-build-debug/src && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /media/sf_CS/git/koi_0x4249/SS20/WR/sc_heat_equation/src/heat.cpp > CMakeFiles/diffusion.dir/heat.cpp.i

src/CMakeFiles/diffusion.dir/heat.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/diffusion.dir/heat.cpp.s"
	cd /media/sf_CS/git/koi_0x4249/SS20/WR/sc_heat_equation/cmake-build-debug/src && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /media/sf_CS/git/koi_0x4249/SS20/WR/sc_heat_equation/src/heat.cpp -o CMakeFiles/diffusion.dir/heat.cpp.s

# Object files for target diffusion
diffusion_OBJECTS = \
"CMakeFiles/diffusion.dir/HeatEquationViewer.cpp.o" \
"CMakeFiles/diffusion.dir/heat.cpp.o"

# External object files for target diffusion
diffusion_EXTERNAL_OBJECTS =

diffusion: src/CMakeFiles/diffusion.dir/HeatEquationViewer.cpp.o
diffusion: src/CMakeFiles/diffusion.dir/heat.cpp.o
diffusion: src/CMakeFiles/diffusion.dir/build.make
diffusion: libglew.a
diffusion: libpmp.a
diffusion: /usr/lib/x86_64-linux-gnu/libGL.so
diffusion: /usr/lib/x86_64-linux-gnu/libGLU.so
diffusion: libglew.a
diffusion: libimgui.a
diffusion: libglfw3.a
diffusion: /usr/lib/x86_64-linux-gnu/librt.so
diffusion: /usr/lib/x86_64-linux-gnu/libm.so
diffusion: /usr/lib/x86_64-linux-gnu/libX11.so
diffusion: src/CMakeFiles/diffusion.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/media/sf_CS/git/koi_0x4249/SS20/WR/sc_heat_equation/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Linking CXX executable ../diffusion"
	cd /media/sf_CS/git/koi_0x4249/SS20/WR/sc_heat_equation/cmake-build-debug/src && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/diffusion.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
src/CMakeFiles/diffusion.dir/build: diffusion

.PHONY : src/CMakeFiles/diffusion.dir/build

src/CMakeFiles/diffusion.dir/clean:
	cd /media/sf_CS/git/koi_0x4249/SS20/WR/sc_heat_equation/cmake-build-debug/src && $(CMAKE_COMMAND) -P CMakeFiles/diffusion.dir/cmake_clean.cmake
.PHONY : src/CMakeFiles/diffusion.dir/clean

src/CMakeFiles/diffusion.dir/depend:
	cd /media/sf_CS/git/koi_0x4249/SS20/WR/sc_heat_equation/cmake-build-debug && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /media/sf_CS/git/koi_0x4249/SS20/WR/sc_heat_equation /media/sf_CS/git/koi_0x4249/SS20/WR/sc_heat_equation/src /media/sf_CS/git/koi_0x4249/SS20/WR/sc_heat_equation/cmake-build-debug /media/sf_CS/git/koi_0x4249/SS20/WR/sc_heat_equation/cmake-build-debug/src /media/sf_CS/git/koi_0x4249/SS20/WR/sc_heat_equation/cmake-build-debug/src/CMakeFiles/diffusion.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : src/CMakeFiles/diffusion.dir/depend

