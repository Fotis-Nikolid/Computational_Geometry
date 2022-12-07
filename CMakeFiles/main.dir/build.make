# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.22

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
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/fotisnikolidais/Desktop/project/Computational_Geometry

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/fotisnikolidais/Desktop/project/Computational_Geometry

# Include any dependencies generated for this target.
include CMakeFiles/main.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/main.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/main.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/main.dir/flags.make

CMakeFiles/main.dir/hull.cc.o: CMakeFiles/main.dir/flags.make
CMakeFiles/main.dir/hull.cc.o: hull.cc
CMakeFiles/main.dir/hull.cc.o: CMakeFiles/main.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/fotisnikolidais/Desktop/project/Computational_Geometry/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/main.dir/hull.cc.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/main.dir/hull.cc.o -MF CMakeFiles/main.dir/hull.cc.o.d -o CMakeFiles/main.dir/hull.cc.o -c /home/fotisnikolidais/Desktop/project/Computational_Geometry/hull.cc

CMakeFiles/main.dir/hull.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/main.dir/hull.cc.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/fotisnikolidais/Desktop/project/Computational_Geometry/hull.cc > CMakeFiles/main.dir/hull.cc.i

CMakeFiles/main.dir/hull.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/main.dir/hull.cc.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/fotisnikolidais/Desktop/project/Computational_Geometry/hull.cc -o CMakeFiles/main.dir/hull.cc.s

CMakeFiles/main.dir/incremental.cpp.o: CMakeFiles/main.dir/flags.make
CMakeFiles/main.dir/incremental.cpp.o: incremental.cpp
CMakeFiles/main.dir/incremental.cpp.o: CMakeFiles/main.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/fotisnikolidais/Desktop/project/Computational_Geometry/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/main.dir/incremental.cpp.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/main.dir/incremental.cpp.o -MF CMakeFiles/main.dir/incremental.cpp.o.d -o CMakeFiles/main.dir/incremental.cpp.o -c /home/fotisnikolidais/Desktop/project/Computational_Geometry/incremental.cpp

CMakeFiles/main.dir/incremental.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/main.dir/incremental.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/fotisnikolidais/Desktop/project/Computational_Geometry/incremental.cpp > CMakeFiles/main.dir/incremental.cpp.i

CMakeFiles/main.dir/incremental.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/main.dir/incremental.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/fotisnikolidais/Desktop/project/Computational_Geometry/incremental.cpp -o CMakeFiles/main.dir/incremental.cpp.s

CMakeFiles/main.dir/kd_tree.cpp.o: CMakeFiles/main.dir/flags.make
CMakeFiles/main.dir/kd_tree.cpp.o: kd_tree.cpp
CMakeFiles/main.dir/kd_tree.cpp.o: CMakeFiles/main.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/fotisnikolidais/Desktop/project/Computational_Geometry/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object CMakeFiles/main.dir/kd_tree.cpp.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/main.dir/kd_tree.cpp.o -MF CMakeFiles/main.dir/kd_tree.cpp.o.d -o CMakeFiles/main.dir/kd_tree.cpp.o -c /home/fotisnikolidais/Desktop/project/Computational_Geometry/kd_tree.cpp

CMakeFiles/main.dir/kd_tree.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/main.dir/kd_tree.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/fotisnikolidais/Desktop/project/Computational_Geometry/kd_tree.cpp > CMakeFiles/main.dir/kd_tree.cpp.i

CMakeFiles/main.dir/kd_tree.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/main.dir/kd_tree.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/fotisnikolidais/Desktop/project/Computational_Geometry/kd_tree.cpp -o CMakeFiles/main.dir/kd_tree.cpp.s

CMakeFiles/main.dir/localsearch.cpp.o: CMakeFiles/main.dir/flags.make
CMakeFiles/main.dir/localsearch.cpp.o: localsearch.cpp
CMakeFiles/main.dir/localsearch.cpp.o: CMakeFiles/main.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/fotisnikolidais/Desktop/project/Computational_Geometry/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object CMakeFiles/main.dir/localsearch.cpp.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/main.dir/localsearch.cpp.o -MF CMakeFiles/main.dir/localsearch.cpp.o.d -o CMakeFiles/main.dir/localsearch.cpp.o -c /home/fotisnikolidais/Desktop/project/Computational_Geometry/localsearch.cpp

CMakeFiles/main.dir/localsearch.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/main.dir/localsearch.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/fotisnikolidais/Desktop/project/Computational_Geometry/localsearch.cpp > CMakeFiles/main.dir/localsearch.cpp.i

CMakeFiles/main.dir/localsearch.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/main.dir/localsearch.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/fotisnikolidais/Desktop/project/Computational_Geometry/localsearch.cpp -o CMakeFiles/main.dir/localsearch.cpp.s

CMakeFiles/main.dir/main.cpp.o: CMakeFiles/main.dir/flags.make
CMakeFiles/main.dir/main.cpp.o: main.cpp
CMakeFiles/main.dir/main.cpp.o: CMakeFiles/main.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/fotisnikolidais/Desktop/project/Computational_Geometry/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building CXX object CMakeFiles/main.dir/main.cpp.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/main.dir/main.cpp.o -MF CMakeFiles/main.dir/main.cpp.o.d -o CMakeFiles/main.dir/main.cpp.o -c /home/fotisnikolidais/Desktop/project/Computational_Geometry/main.cpp

CMakeFiles/main.dir/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/main.dir/main.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/fotisnikolidais/Desktop/project/Computational_Geometry/main.cpp > CMakeFiles/main.dir/main.cpp.i

CMakeFiles/main.dir/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/main.dir/main.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/fotisnikolidais/Desktop/project/Computational_Geometry/main.cpp -o CMakeFiles/main.dir/main.cpp.s

CMakeFiles/main.dir/polygon.cc.o: CMakeFiles/main.dir/flags.make
CMakeFiles/main.dir/polygon.cc.o: polygon.cc
CMakeFiles/main.dir/polygon.cc.o: CMakeFiles/main.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/fotisnikolidais/Desktop/project/Computational_Geometry/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Building CXX object CMakeFiles/main.dir/polygon.cc.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/main.dir/polygon.cc.o -MF CMakeFiles/main.dir/polygon.cc.o.d -o CMakeFiles/main.dir/polygon.cc.o -c /home/fotisnikolidais/Desktop/project/Computational_Geometry/polygon.cc

CMakeFiles/main.dir/polygon.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/main.dir/polygon.cc.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/fotisnikolidais/Desktop/project/Computational_Geometry/polygon.cc > CMakeFiles/main.dir/polygon.cc.i

CMakeFiles/main.dir/polygon.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/main.dir/polygon.cc.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/fotisnikolidais/Desktop/project/Computational_Geometry/polygon.cc -o CMakeFiles/main.dir/polygon.cc.s

CMakeFiles/main.dir/simulated_annealing.cpp.o: CMakeFiles/main.dir/flags.make
CMakeFiles/main.dir/simulated_annealing.cpp.o: simulated_annealing.cpp
CMakeFiles/main.dir/simulated_annealing.cpp.o: CMakeFiles/main.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/fotisnikolidais/Desktop/project/Computational_Geometry/CMakeFiles --progress-num=$(CMAKE_PROGRESS_7) "Building CXX object CMakeFiles/main.dir/simulated_annealing.cpp.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/main.dir/simulated_annealing.cpp.o -MF CMakeFiles/main.dir/simulated_annealing.cpp.o.d -o CMakeFiles/main.dir/simulated_annealing.cpp.o -c /home/fotisnikolidais/Desktop/project/Computational_Geometry/simulated_annealing.cpp

CMakeFiles/main.dir/simulated_annealing.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/main.dir/simulated_annealing.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/fotisnikolidais/Desktop/project/Computational_Geometry/simulated_annealing.cpp > CMakeFiles/main.dir/simulated_annealing.cpp.i

CMakeFiles/main.dir/simulated_annealing.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/main.dir/simulated_annealing.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/fotisnikolidais/Desktop/project/Computational_Geometry/simulated_annealing.cpp -o CMakeFiles/main.dir/simulated_annealing.cpp.s

CMakeFiles/main.dir/test.cpp.o: CMakeFiles/main.dir/flags.make
CMakeFiles/main.dir/test.cpp.o: test.cpp
CMakeFiles/main.dir/test.cpp.o: CMakeFiles/main.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/fotisnikolidais/Desktop/project/Computational_Geometry/CMakeFiles --progress-num=$(CMAKE_PROGRESS_8) "Building CXX object CMakeFiles/main.dir/test.cpp.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/main.dir/test.cpp.o -MF CMakeFiles/main.dir/test.cpp.o.d -o CMakeFiles/main.dir/test.cpp.o -c /home/fotisnikolidais/Desktop/project/Computational_Geometry/test.cpp

CMakeFiles/main.dir/test.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/main.dir/test.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/fotisnikolidais/Desktop/project/Computational_Geometry/test.cpp > CMakeFiles/main.dir/test.cpp.i

CMakeFiles/main.dir/test.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/main.dir/test.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/fotisnikolidais/Desktop/project/Computational_Geometry/test.cpp -o CMakeFiles/main.dir/test.cpp.s

# Object files for target main
main_OBJECTS = \
"CMakeFiles/main.dir/hull.cc.o" \
"CMakeFiles/main.dir/incremental.cpp.o" \
"CMakeFiles/main.dir/kd_tree.cpp.o" \
"CMakeFiles/main.dir/localsearch.cpp.o" \
"CMakeFiles/main.dir/main.cpp.o" \
"CMakeFiles/main.dir/polygon.cc.o" \
"CMakeFiles/main.dir/simulated_annealing.cpp.o" \
"CMakeFiles/main.dir/test.cpp.o"

# External object files for target main
main_EXTERNAL_OBJECTS =

main: CMakeFiles/main.dir/hull.cc.o
main: CMakeFiles/main.dir/incremental.cpp.o
main: CMakeFiles/main.dir/kd_tree.cpp.o
main: CMakeFiles/main.dir/localsearch.cpp.o
main: CMakeFiles/main.dir/main.cpp.o
main: CMakeFiles/main.dir/polygon.cc.o
main: CMakeFiles/main.dir/simulated_annealing.cpp.o
main: CMakeFiles/main.dir/test.cpp.o
main: CMakeFiles/main.dir/build.make
main: /usr/lib/x86_64-linux-gnu/libgmpxx.so
main: /usr/lib/x86_64-linux-gnu/libmpfr.so
main: /usr/lib/x86_64-linux-gnu/libgmp.so
main: CMakeFiles/main.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/fotisnikolidais/Desktop/project/Computational_Geometry/CMakeFiles --progress-num=$(CMAKE_PROGRESS_9) "Linking CXX executable main"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/main.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/main.dir/build: main
.PHONY : CMakeFiles/main.dir/build

CMakeFiles/main.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/main.dir/cmake_clean.cmake
.PHONY : CMakeFiles/main.dir/clean

CMakeFiles/main.dir/depend:
	cd /home/fotisnikolidais/Desktop/project/Computational_Geometry && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/fotisnikolidais/Desktop/project/Computational_Geometry /home/fotisnikolidais/Desktop/project/Computational_Geometry /home/fotisnikolidais/Desktop/project/Computational_Geometry /home/fotisnikolidais/Desktop/project/Computational_Geometry /home/fotisnikolidais/Desktop/project/Computational_Geometry/CMakeFiles/main.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/main.dir/depend

