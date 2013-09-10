# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 2.8

# Default target executed when no arguments are given to make.
default_target: all
.PHONY : default_target

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canoncical targets will work.
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

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/dimmockm/Private/workspace/geant/dxfm

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/dimmockm/Private/workspace/geant/dxfm-build

#=============================================================================
# Targets provided globally by CMake.

# Special rule for the target edit_cache
edit_cache:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Running interactive CMake command-line interface..."
	/usr/bin/cmake -i .
.PHONY : edit_cache

# Special rule for the target edit_cache
edit_cache/fast: edit_cache
.PHONY : edit_cache/fast

# Special rule for the target install
install: preinstall
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Install the project..."
	/usr/bin/cmake -P cmake_install.cmake
.PHONY : install

# Special rule for the target install
install/fast: preinstall/fast
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Install the project..."
	/usr/bin/cmake -P cmake_install.cmake
.PHONY : install/fast

# Special rule for the target install/local
install/local: preinstall
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Installing only the local directory..."
	/usr/bin/cmake -DCMAKE_INSTALL_LOCAL_ONLY=1 -P cmake_install.cmake
.PHONY : install/local

# Special rule for the target install/local
install/local/fast: install/local
.PHONY : install/local/fast

# Special rule for the target install/strip
install/strip: preinstall
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Installing the project stripped..."
	/usr/bin/cmake -DCMAKE_INSTALL_DO_STRIP=1 -P cmake_install.cmake
.PHONY : install/strip

# Special rule for the target install/strip
install/strip/fast: install/strip
.PHONY : install/strip/fast

# Special rule for the target list_install_components
list_install_components:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Available install components are: \"Unspecified\""
.PHONY : list_install_components

# Special rule for the target list_install_components
list_install_components/fast: list_install_components
.PHONY : list_install_components/fast

# Special rule for the target rebuild_cache
rebuild_cache:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Running CMake to regenerate build system..."
	/usr/bin/cmake -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR)
.PHONY : rebuild_cache

# Special rule for the target rebuild_cache
rebuild_cache/fast: rebuild_cache
.PHONY : rebuild_cache/fast

# The main all target
all: cmake_check_build_system
	$(CMAKE_COMMAND) -E cmake_progress_start /home/dimmockm/Private/workspace/geant/dxfm-build/CMakeFiles /home/dimmockm/Private/workspace/geant/dxfm-build/CMakeFiles/progress.marks
	$(MAKE) -f CMakeFiles/Makefile2 all
	$(CMAKE_COMMAND) -E cmake_progress_start /home/dimmockm/Private/workspace/geant/dxfm-build/CMakeFiles 0
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
	$(CMAKE_COMMAND) -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR) --check-build-system CMakeFiles/Makefile.cmake 1
.PHONY : depend

#=============================================================================
# Target rules for targets named dxfm

# Build rule for target.
dxfm: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 dxfm
.PHONY : dxfm

# fast build rule for target.
dxfm/fast:
	$(MAKE) -f CMakeFiles/dxfm.dir/build.make CMakeFiles/dxfm.dir/build
.PHONY : dxfm/fast

#=============================================================================
# Target rules for targets named fluo

# Build rule for target.
fluo: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 fluo
.PHONY : fluo

# fast build rule for target.
fluo/fast:
	$(MAKE) -f CMakeFiles/fluo.dir/build.make CMakeFiles/fluo.dir/build
.PHONY : fluo/fast

fluoro.o: fluoro.cc.o
.PHONY : fluoro.o

# target to build an object file
fluoro.cc.o:
	$(MAKE) -f CMakeFiles/fluo.dir/build.make CMakeFiles/fluo.dir/fluoro.cc.o
.PHONY : fluoro.cc.o

fluoro.i: fluoro.cc.i
.PHONY : fluoro.i

# target to preprocess a source file
fluoro.cc.i:
	$(MAKE) -f CMakeFiles/fluo.dir/build.make CMakeFiles/fluo.dir/fluoro.cc.i
.PHONY : fluoro.cc.i

fluoro.s: fluoro.cc.s
.PHONY : fluoro.s

# target to generate assembly for a file
fluoro.cc.s:
	$(MAKE) -f CMakeFiles/fluo.dir/build.make CMakeFiles/fluo.dir/fluoro.cc.s
.PHONY : fluoro.cc.s

src/DetectorConstruction.o: src/DetectorConstruction.cc.o
.PHONY : src/DetectorConstruction.o

# target to build an object file
src/DetectorConstruction.cc.o:
	$(MAKE) -f CMakeFiles/fluo.dir/build.make CMakeFiles/fluo.dir/src/DetectorConstruction.cc.o
.PHONY : src/DetectorConstruction.cc.o

src/DetectorConstruction.i: src/DetectorConstruction.cc.i
.PHONY : src/DetectorConstruction.i

# target to preprocess a source file
src/DetectorConstruction.cc.i:
	$(MAKE) -f CMakeFiles/fluo.dir/build.make CMakeFiles/fluo.dir/src/DetectorConstruction.cc.i
.PHONY : src/DetectorConstruction.cc.i

src/DetectorConstruction.s: src/DetectorConstruction.cc.s
.PHONY : src/DetectorConstruction.s

# target to generate assembly for a file
src/DetectorConstruction.cc.s:
	$(MAKE) -f CMakeFiles/fluo.dir/build.make CMakeFiles/fluo.dir/src/DetectorConstruction.cc.s
.PHONY : src/DetectorConstruction.cc.s

src/DetectorHits.o: src/DetectorHits.cc.o
.PHONY : src/DetectorHits.o

# target to build an object file
src/DetectorHits.cc.o:
	$(MAKE) -f CMakeFiles/fluo.dir/build.make CMakeFiles/fluo.dir/src/DetectorHits.cc.o
.PHONY : src/DetectorHits.cc.o

src/DetectorHits.i: src/DetectorHits.cc.i
.PHONY : src/DetectorHits.i

# target to preprocess a source file
src/DetectorHits.cc.i:
	$(MAKE) -f CMakeFiles/fluo.dir/build.make CMakeFiles/fluo.dir/src/DetectorHits.cc.i
.PHONY : src/DetectorHits.cc.i

src/DetectorHits.s: src/DetectorHits.cc.s
.PHONY : src/DetectorHits.s

# target to generate assembly for a file
src/DetectorHits.cc.s:
	$(MAKE) -f CMakeFiles/fluo.dir/build.make CMakeFiles/fluo.dir/src/DetectorHits.cc.s
.PHONY : src/DetectorHits.cc.s

src/DetectorMessenger.o: src/DetectorMessenger.cc.o
.PHONY : src/DetectorMessenger.o

# target to build an object file
src/DetectorMessenger.cc.o:
	$(MAKE) -f CMakeFiles/fluo.dir/build.make CMakeFiles/fluo.dir/src/DetectorMessenger.cc.o
.PHONY : src/DetectorMessenger.cc.o

src/DetectorMessenger.i: src/DetectorMessenger.cc.i
.PHONY : src/DetectorMessenger.i

# target to preprocess a source file
src/DetectorMessenger.cc.i:
	$(MAKE) -f CMakeFiles/fluo.dir/build.make CMakeFiles/fluo.dir/src/DetectorMessenger.cc.i
.PHONY : src/DetectorMessenger.cc.i

src/DetectorMessenger.s: src/DetectorMessenger.cc.s
.PHONY : src/DetectorMessenger.s

# target to generate assembly for a file
src/DetectorMessenger.cc.s:
	$(MAKE) -f CMakeFiles/fluo.dir/build.make CMakeFiles/fluo.dir/src/DetectorMessenger.cc.s
.PHONY : src/DetectorMessenger.cc.s

src/EventAction.o: src/EventAction.cc.o
.PHONY : src/EventAction.o

# target to build an object file
src/EventAction.cc.o:
	$(MAKE) -f CMakeFiles/fluo.dir/build.make CMakeFiles/fluo.dir/src/EventAction.cc.o
.PHONY : src/EventAction.cc.o

src/EventAction.i: src/EventAction.cc.i
.PHONY : src/EventAction.i

# target to preprocess a source file
src/EventAction.cc.i:
	$(MAKE) -f CMakeFiles/fluo.dir/build.make CMakeFiles/fluo.dir/src/EventAction.cc.i
.PHONY : src/EventAction.cc.i

src/EventAction.s: src/EventAction.cc.s
.PHONY : src/EventAction.s

# target to generate assembly for a file
src/EventAction.cc.s:
	$(MAKE) -f CMakeFiles/fluo.dir/build.make CMakeFiles/fluo.dir/src/EventAction.cc.s
.PHONY : src/EventAction.cc.s

src/EventMessenger.o: src/EventMessenger.cc.o
.PHONY : src/EventMessenger.o

# target to build an object file
src/EventMessenger.cc.o:
	$(MAKE) -f CMakeFiles/fluo.dir/build.make CMakeFiles/fluo.dir/src/EventMessenger.cc.o
.PHONY : src/EventMessenger.cc.o

src/EventMessenger.i: src/EventMessenger.cc.i
.PHONY : src/EventMessenger.i

# target to preprocess a source file
src/EventMessenger.cc.i:
	$(MAKE) -f CMakeFiles/fluo.dir/build.make CMakeFiles/fluo.dir/src/EventMessenger.cc.i
.PHONY : src/EventMessenger.cc.i

src/EventMessenger.s: src/EventMessenger.cc.s
.PHONY : src/EventMessenger.s

# target to generate assembly for a file
src/EventMessenger.cc.s:
	$(MAKE) -f CMakeFiles/fluo.dir/build.make CMakeFiles/fluo.dir/src/EventMessenger.cc.s
.PHONY : src/EventMessenger.cc.s

src/ExpPhantomHit.o: src/ExpPhantomHit.cc.o
.PHONY : src/ExpPhantomHit.o

# target to build an object file
src/ExpPhantomHit.cc.o:
	$(MAKE) -f CMakeFiles/fluo.dir/build.make CMakeFiles/fluo.dir/src/ExpPhantomHit.cc.o
.PHONY : src/ExpPhantomHit.cc.o

src/ExpPhantomHit.i: src/ExpPhantomHit.cc.i
.PHONY : src/ExpPhantomHit.i

# target to preprocess a source file
src/ExpPhantomHit.cc.i:
	$(MAKE) -f CMakeFiles/fluo.dir/build.make CMakeFiles/fluo.dir/src/ExpPhantomHit.cc.i
.PHONY : src/ExpPhantomHit.cc.i

src/ExpPhantomHit.s: src/ExpPhantomHit.cc.s
.PHONY : src/ExpPhantomHit.s

# target to generate assembly for a file
src/ExpPhantomHit.cc.s:
	$(MAKE) -f CMakeFiles/fluo.dir/build.make CMakeFiles/fluo.dir/src/ExpPhantomHit.cc.s
.PHONY : src/ExpPhantomHit.cc.s

src/ExpPhantomSD.o: src/ExpPhantomSD.cc.o
.PHONY : src/ExpPhantomSD.o

# target to build an object file
src/ExpPhantomSD.cc.o:
	$(MAKE) -f CMakeFiles/fluo.dir/build.make CMakeFiles/fluo.dir/src/ExpPhantomSD.cc.o
.PHONY : src/ExpPhantomSD.cc.o

src/ExpPhantomSD.i: src/ExpPhantomSD.cc.i
.PHONY : src/ExpPhantomSD.i

# target to preprocess a source file
src/ExpPhantomSD.cc.i:
	$(MAKE) -f CMakeFiles/fluo.dir/build.make CMakeFiles/fluo.dir/src/ExpPhantomSD.cc.i
.PHONY : src/ExpPhantomSD.cc.i

src/ExpPhantomSD.s: src/ExpPhantomSD.cc.s
.PHONY : src/ExpPhantomSD.s

# target to generate assembly for a file
src/ExpPhantomSD.cc.s:
	$(MAKE) -f CMakeFiles/fluo.dir/build.make CMakeFiles/fluo.dir/src/ExpPhantomSD.cc.s
.PHONY : src/ExpPhantomSD.cc.s

src/PhysicsList.o: src/PhysicsList.cc.o
.PHONY : src/PhysicsList.o

# target to build an object file
src/PhysicsList.cc.o:
	$(MAKE) -f CMakeFiles/fluo.dir/build.make CMakeFiles/fluo.dir/src/PhysicsList.cc.o
.PHONY : src/PhysicsList.cc.o

src/PhysicsList.i: src/PhysicsList.cc.i
.PHONY : src/PhysicsList.i

# target to preprocess a source file
src/PhysicsList.cc.i:
	$(MAKE) -f CMakeFiles/fluo.dir/build.make CMakeFiles/fluo.dir/src/PhysicsList.cc.i
.PHONY : src/PhysicsList.cc.i

src/PhysicsList.s: src/PhysicsList.cc.s
.PHONY : src/PhysicsList.s

# target to generate assembly for a file
src/PhysicsList.cc.s:
	$(MAKE) -f CMakeFiles/fluo.dir/build.make CMakeFiles/fluo.dir/src/PhysicsList.cc.s
.PHONY : src/PhysicsList.cc.s

src/PrimaryGeneratorAction.o: src/PrimaryGeneratorAction.cc.o
.PHONY : src/PrimaryGeneratorAction.o

# target to build an object file
src/PrimaryGeneratorAction.cc.o:
	$(MAKE) -f CMakeFiles/fluo.dir/build.make CMakeFiles/fluo.dir/src/PrimaryGeneratorAction.cc.o
.PHONY : src/PrimaryGeneratorAction.cc.o

src/PrimaryGeneratorAction.i: src/PrimaryGeneratorAction.cc.i
.PHONY : src/PrimaryGeneratorAction.i

# target to preprocess a source file
src/PrimaryGeneratorAction.cc.i:
	$(MAKE) -f CMakeFiles/fluo.dir/build.make CMakeFiles/fluo.dir/src/PrimaryGeneratorAction.cc.i
.PHONY : src/PrimaryGeneratorAction.cc.i

src/PrimaryGeneratorAction.s: src/PrimaryGeneratorAction.cc.s
.PHONY : src/PrimaryGeneratorAction.s

# target to generate assembly for a file
src/PrimaryGeneratorAction.cc.s:
	$(MAKE) -f CMakeFiles/fluo.dir/build.make CMakeFiles/fluo.dir/src/PrimaryGeneratorAction.cc.s
.PHONY : src/PrimaryGeneratorAction.cc.s

src/PrimaryGeneratorMessenger.o: src/PrimaryGeneratorMessenger.cc.o
.PHONY : src/PrimaryGeneratorMessenger.o

# target to build an object file
src/PrimaryGeneratorMessenger.cc.o:
	$(MAKE) -f CMakeFiles/fluo.dir/build.make CMakeFiles/fluo.dir/src/PrimaryGeneratorMessenger.cc.o
.PHONY : src/PrimaryGeneratorMessenger.cc.o

src/PrimaryGeneratorMessenger.i: src/PrimaryGeneratorMessenger.cc.i
.PHONY : src/PrimaryGeneratorMessenger.i

# target to preprocess a source file
src/PrimaryGeneratorMessenger.cc.i:
	$(MAKE) -f CMakeFiles/fluo.dir/build.make CMakeFiles/fluo.dir/src/PrimaryGeneratorMessenger.cc.i
.PHONY : src/PrimaryGeneratorMessenger.cc.i

src/PrimaryGeneratorMessenger.s: src/PrimaryGeneratorMessenger.cc.s
.PHONY : src/PrimaryGeneratorMessenger.s

# target to generate assembly for a file
src/PrimaryGeneratorMessenger.cc.s:
	$(MAKE) -f CMakeFiles/fluo.dir/build.make CMakeFiles/fluo.dir/src/PrimaryGeneratorMessenger.cc.s
.PHONY : src/PrimaryGeneratorMessenger.cc.s

src/RunAction.o: src/RunAction.cc.o
.PHONY : src/RunAction.o

# target to build an object file
src/RunAction.cc.o:
	$(MAKE) -f CMakeFiles/fluo.dir/build.make CMakeFiles/fluo.dir/src/RunAction.cc.o
.PHONY : src/RunAction.cc.o

src/RunAction.i: src/RunAction.cc.i
.PHONY : src/RunAction.i

# target to preprocess a source file
src/RunAction.cc.i:
	$(MAKE) -f CMakeFiles/fluo.dir/build.make CMakeFiles/fluo.dir/src/RunAction.cc.i
.PHONY : src/RunAction.cc.i

src/RunAction.s: src/RunAction.cc.s
.PHONY : src/RunAction.s

# target to generate assembly for a file
src/RunAction.cc.s:
	$(MAKE) -f CMakeFiles/fluo.dir/build.make CMakeFiles/fluo.dir/src/RunAction.cc.s
.PHONY : src/RunAction.cc.s

src/RunMessenger.o: src/RunMessenger.cc.o
.PHONY : src/RunMessenger.o

# target to build an object file
src/RunMessenger.cc.o:
	$(MAKE) -f CMakeFiles/fluo.dir/build.make CMakeFiles/fluo.dir/src/RunMessenger.cc.o
.PHONY : src/RunMessenger.cc.o

src/RunMessenger.i: src/RunMessenger.cc.i
.PHONY : src/RunMessenger.i

# target to preprocess a source file
src/RunMessenger.cc.i:
	$(MAKE) -f CMakeFiles/fluo.dir/build.make CMakeFiles/fluo.dir/src/RunMessenger.cc.i
.PHONY : src/RunMessenger.cc.i

src/RunMessenger.s: src/RunMessenger.cc.s
.PHONY : src/RunMessenger.s

# target to generate assembly for a file
src/RunMessenger.cc.s:
	$(MAKE) -f CMakeFiles/fluo.dir/build.make CMakeFiles/fluo.dir/src/RunMessenger.cc.s
.PHONY : src/RunMessenger.cc.s

src/SensitiveDet.o: src/SensitiveDet.cc.o
.PHONY : src/SensitiveDet.o

# target to build an object file
src/SensitiveDet.cc.o:
	$(MAKE) -f CMakeFiles/fluo.dir/build.make CMakeFiles/fluo.dir/src/SensitiveDet.cc.o
.PHONY : src/SensitiveDet.cc.o

src/SensitiveDet.i: src/SensitiveDet.cc.i
.PHONY : src/SensitiveDet.i

# target to preprocess a source file
src/SensitiveDet.cc.i:
	$(MAKE) -f CMakeFiles/fluo.dir/build.make CMakeFiles/fluo.dir/src/SensitiveDet.cc.i
.PHONY : src/SensitiveDet.cc.i

src/SensitiveDet.s: src/SensitiveDet.cc.s
.PHONY : src/SensitiveDet.s

# target to generate assembly for a file
src/SensitiveDet.cc.s:
	$(MAKE) -f CMakeFiles/fluo.dir/build.make CMakeFiles/fluo.dir/src/SensitiveDet.cc.s
.PHONY : src/SensitiveDet.cc.s

src/VisManager.o: src/VisManager.cc.o
.PHONY : src/VisManager.o

# target to build an object file
src/VisManager.cc.o:
	$(MAKE) -f CMakeFiles/fluo.dir/build.make CMakeFiles/fluo.dir/src/VisManager.cc.o
.PHONY : src/VisManager.cc.o

src/VisManager.i: src/VisManager.cc.i
.PHONY : src/VisManager.i

# target to preprocess a source file
src/VisManager.cc.i:
	$(MAKE) -f CMakeFiles/fluo.dir/build.make CMakeFiles/fluo.dir/src/VisManager.cc.i
.PHONY : src/VisManager.cc.i

src/VisManager.s: src/VisManager.cc.s
.PHONY : src/VisManager.s

# target to generate assembly for a file
src/VisManager.cc.s:
	$(MAKE) -f CMakeFiles/fluo.dir/build.make CMakeFiles/fluo.dir/src/VisManager.cc.s
.PHONY : src/VisManager.cc.s

# Help Target
help:
	@echo "The following are some of the valid targets for this Makefile:"
	@echo "... all (the default if no target is provided)"
	@echo "... clean"
	@echo "... depend"
	@echo "... dxfm"
	@echo "... edit_cache"
	@echo "... fluo"
	@echo "... install"
	@echo "... install/local"
	@echo "... install/strip"
	@echo "... list_install_components"
	@echo "... rebuild_cache"
	@echo "... fluoro.o"
	@echo "... fluoro.i"
	@echo "... fluoro.s"
	@echo "... src/DetectorConstruction.o"
	@echo "... src/DetectorConstruction.i"
	@echo "... src/DetectorConstruction.s"
	@echo "... src/DetectorHits.o"
	@echo "... src/DetectorHits.i"
	@echo "... src/DetectorHits.s"
	@echo "... src/DetectorMessenger.o"
	@echo "... src/DetectorMessenger.i"
	@echo "... src/DetectorMessenger.s"
	@echo "... src/EventAction.o"
	@echo "... src/EventAction.i"
	@echo "... src/EventAction.s"
	@echo "... src/EventMessenger.o"
	@echo "... src/EventMessenger.i"
	@echo "... src/EventMessenger.s"
	@echo "... src/ExpPhantomHit.o"
	@echo "... src/ExpPhantomHit.i"
	@echo "... src/ExpPhantomHit.s"
	@echo "... src/ExpPhantomSD.o"
	@echo "... src/ExpPhantomSD.i"
	@echo "... src/ExpPhantomSD.s"
	@echo "... src/PhysicsList.o"
	@echo "... src/PhysicsList.i"
	@echo "... src/PhysicsList.s"
	@echo "... src/PrimaryGeneratorAction.o"
	@echo "... src/PrimaryGeneratorAction.i"
	@echo "... src/PrimaryGeneratorAction.s"
	@echo "... src/PrimaryGeneratorMessenger.o"
	@echo "... src/PrimaryGeneratorMessenger.i"
	@echo "... src/PrimaryGeneratorMessenger.s"
	@echo "... src/RunAction.o"
	@echo "... src/RunAction.i"
	@echo "... src/RunAction.s"
	@echo "... src/RunMessenger.o"
	@echo "... src/RunMessenger.i"
	@echo "... src/RunMessenger.s"
	@echo "... src/SensitiveDet.o"
	@echo "... src/SensitiveDet.i"
	@echo "... src/SensitiveDet.s"
	@echo "... src/VisManager.o"
	@echo "... src/VisManager.i"
	@echo "... src/VisManager.s"
.PHONY : help



#=============================================================================
# Special targets to cleanup operation of make.

# Special rule to run CMake to check the build system integrity.
# No rule that depends on this can have commands that come from listfiles
# because they might be regenerated.
cmake_check_build_system:
	$(CMAKE_COMMAND) -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR) --check-build-system CMakeFiles/Makefile.cmake 0
.PHONY : cmake_check_build_system

