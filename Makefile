# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.16

# Default target executed when no arguments are given to make.
default_target: all

.PHONY : default_target

# Allow only one "make -f Makefile2" at a time, but pass parallelism.
.NOTPARALLEL:


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
CMAKE_SOURCE_DIR = /home/users/sdi1900132/p1final/Computational_Geometry

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/users/sdi1900132/p1final/Computational_Geometry

#=============================================================================
# Targets provided globally by CMake.

# Special rule for the target rebuild_cache
rebuild_cache:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Running CMake to regenerate build system..."
	/usr/bin/cmake -S$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR)
.PHONY : rebuild_cache

# Special rule for the target rebuild_cache
rebuild_cache/fast: rebuild_cache

.PHONY : rebuild_cache/fast

# Special rule for the target edit_cache
edit_cache:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "No interactive CMake dialog available..."
	/usr/bin/cmake -E echo No\ interactive\ CMake\ dialog\ available.
.PHONY : edit_cache

# Special rule for the target edit_cache
edit_cache/fast: edit_cache

.PHONY : edit_cache/fast

# The main all target
all: cmake_check_build_system
	$(CMAKE_COMMAND) -E cmake_progress_start /home/users/sdi1900132/p1final/Computational_Geometry/CMakeFiles /home/users/sdi1900132/p1final/Computational_Geometry/CMakeFiles/progress.marks
	$(MAKE) -f CMakeFiles/Makefile2 all
	$(CMAKE_COMMAND) -E cmake_progress_start /home/users/sdi1900132/p1final/Computational_Geometry/CMakeFiles 0
.PHONY : all

# The main clean target
clean:
	$(MAKE) -f CMakeFiles/Makefile2 clean
.PHONY : clean

# The main clean target
clean/fast: clean

.PHONY : clean/fast

# Prepare targets for installation.
preinstall: all
	$(MAKE) -f CMakeFiles/Makefile2 preinstall
.PHONY : preinstall

# Prepare targets for installation.
preinstall/fast:
	$(MAKE) -f CMakeFiles/Makefile2 preinstall
.PHONY : preinstall/fast

# clear depends
depend:
	$(CMAKE_COMMAND) -S$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR) --check-build-system CMakeFiles/Makefile.cmake 1
.PHONY : depend

#=============================================================================
# Target rules for targets named to_polygon

# Build rule for target.
to_polygon: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 to_polygon
.PHONY : to_polygon

# fast build rule for target.
to_polygon/fast:
	$(MAKE) -f CMakeFiles/to_polygon.dir/build.make CMakeFiles/to_polygon.dir/build
.PHONY : to_polygon/fast

hull.o: hull.cc.o

.PHONY : hull.o

# target to build an object file
hull.cc.o:
	$(MAKE) -f CMakeFiles/to_polygon.dir/build.make CMakeFiles/to_polygon.dir/hull.cc.o
.PHONY : hull.cc.o

hull.i: hull.cc.i

.PHONY : hull.i

# target to preprocess a source file
hull.cc.i:
	$(MAKE) -f CMakeFiles/to_polygon.dir/build.make CMakeFiles/to_polygon.dir/hull.cc.i
.PHONY : hull.cc.i

hull.s: hull.cc.s

.PHONY : hull.s

# target to generate assembly for a file
hull.cc.s:
	$(MAKE) -f CMakeFiles/to_polygon.dir/build.make CMakeFiles/to_polygon.dir/hull.cc.s
.PHONY : hull.cc.s

incremental.o: incremental.cpp.o

.PHONY : incremental.o

# target to build an object file
incremental.cpp.o:
	$(MAKE) -f CMakeFiles/to_polygon.dir/build.make CMakeFiles/to_polygon.dir/incremental.cpp.o
.PHONY : incremental.cpp.o

incremental.i: incremental.cpp.i

.PHONY : incremental.i

# target to preprocess a source file
incremental.cpp.i:
	$(MAKE) -f CMakeFiles/to_polygon.dir/build.make CMakeFiles/to_polygon.dir/incremental.cpp.i
.PHONY : incremental.cpp.i

incremental.s: incremental.cpp.s

.PHONY : incremental.s

# target to generate assembly for a file
incremental.cpp.s:
	$(MAKE) -f CMakeFiles/to_polygon.dir/build.make CMakeFiles/to_polygon.dir/incremental.cpp.s
.PHONY : incremental.cpp.s

main.o: main.cpp.o

.PHONY : main.o

# target to build an object file
main.cpp.o:
	$(MAKE) -f CMakeFiles/to_polygon.dir/build.make CMakeFiles/to_polygon.dir/main.cpp.o
.PHONY : main.cpp.o

main.i: main.cpp.i

.PHONY : main.i

# target to preprocess a source file
main.cpp.i:
	$(MAKE) -f CMakeFiles/to_polygon.dir/build.make CMakeFiles/to_polygon.dir/main.cpp.i
.PHONY : main.cpp.i

main.s: main.cpp.s

.PHONY : main.s

# target to generate assembly for a file
main.cpp.s:
	$(MAKE) -f CMakeFiles/to_polygon.dir/build.make CMakeFiles/to_polygon.dir/main.cpp.s
.PHONY : main.cpp.s

polygon.o: polygon.cc.o

.PHONY : polygon.o

# target to build an object file
polygon.cc.o:
	$(MAKE) -f CMakeFiles/to_polygon.dir/build.make CMakeFiles/to_polygon.dir/polygon.cc.o
.PHONY : polygon.cc.o

polygon.i: polygon.cc.i

.PHONY : polygon.i

# target to preprocess a source file
polygon.cc.i:
	$(MAKE) -f CMakeFiles/to_polygon.dir/build.make CMakeFiles/to_polygon.dir/polygon.cc.i
.PHONY : polygon.cc.i

polygon.s: polygon.cc.s

.PHONY : polygon.s

# target to generate assembly for a file
polygon.cc.s:
	$(MAKE) -f CMakeFiles/to_polygon.dir/build.make CMakeFiles/to_polygon.dir/polygon.cc.s
.PHONY : polygon.cc.s

test.o: test.cpp.o

.PHONY : test.o

# target to build an object file
test.cpp.o:
	$(MAKE) -f CMakeFiles/to_polygon.dir/build.make CMakeFiles/to_polygon.dir/test.cpp.o
.PHONY : test.cpp.o

test.i: test.cpp.i

.PHONY : test.i

# target to preprocess a source file
test.cpp.i:
	$(MAKE) -f CMakeFiles/to_polygon.dir/build.make CMakeFiles/to_polygon.dir/test.cpp.i
.PHONY : test.cpp.i

test.s: test.cpp.s

.PHONY : test.s

# target to generate assembly for a file
test.cpp.s:
	$(MAKE) -f CMakeFiles/to_polygon.dir/build.make CMakeFiles/to_polygon.dir/test.cpp.s
.PHONY : test.cpp.s

# Help Target
help:
	@echo "The following are some of the valid targets for this Makefile:"
	@echo "... all (the default if no target is provided)"
	@echo "... clean"
	@echo "... depend"
	@echo "... rebuild_cache"
	@echo "... edit_cache"
	@echo "... to_polygon"
	@echo "... hull.o"
	@echo "... hull.i"
	@echo "... hull.s"
	@echo "... incremental.o"
	@echo "... incremental.i"
	@echo "... incremental.s"
	@echo "... main.o"
	@echo "... main.i"
	@echo "... main.s"
	@echo "... polygon.o"
	@echo "... polygon.i"
	@echo "... polygon.s"
	@echo "... test.o"
	@echo "... test.i"
	@echo "... test.s"
.PHONY : help



#=============================================================================
# Special targets to cleanup operation of make.

# Special rule to run CMake to check the build system integrity.
# No rule that depends on this can have commands that come from listfiles
# because they might be regenerated.
cmake_check_build_system:
	$(CMAKE_COMMAND) -S$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR) --check-build-system CMakeFiles/Makefile.cmake 0
.PHONY : cmake_check_build_system

