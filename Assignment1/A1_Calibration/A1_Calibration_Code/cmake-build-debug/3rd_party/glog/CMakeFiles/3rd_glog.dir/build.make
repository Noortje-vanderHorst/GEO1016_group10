# CMAKE generated file: DO NOT EDIT!
# Generated by "NMake Makefiles" Generator, CMake Version 3.17

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


.SUFFIXES: .hpux_make_needs_suffix_list


# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

!IF "$(OS)" == "Windows_NT"
NULL=
!ELSE
NULL=nul
!ENDIF
SHELL = cmd.exe

# The CMake executable.
CMAKE_COMMAND = "C:\Program Files\JetBrains\CLion 2020.3\bin\cmake\win\bin\cmake.exe"

# The command to remove a file.
RM = "C:\Program Files\JetBrains\CLion 2020.3\bin\cmake\win\bin\cmake.exe" -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = C:\Users\noort\Documents\1\TU\year_4\GEO1016_3D_computer_vision\Assignments\GEO1016_group10\Assignment1\A1_Calibration\A1_Calibration_Code

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = C:\Users\noort\Documents\1\TU\year_4\GEO1016_3D_computer_vision\Assignments\GEO1016_group10\Assignment1\A1_Calibration\A1_Calibration_Code\cmake-build-debug

# Include any dependencies generated for this target.
include 3rd_party\glog\CMakeFiles\3rd_glog.dir\depend.make

# Include the progress variables for this target.
include 3rd_party\glog\CMakeFiles\3rd_glog.dir\progress.make

# Include the compile flags for this target's objects.
include 3rd_party\glog\CMakeFiles\3rd_glog.dir\flags.make

3rd_party\glog\CMakeFiles\3rd_glog.dir\logging.cc.obj: 3rd_party\glog\CMakeFiles\3rd_glog.dir\flags.make
3rd_party\glog\CMakeFiles\3rd_glog.dir\logging.cc.obj: ..\3rd_party\glog\logging.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=C:\Users\noort\Documents\1\TU\year_4\GEO1016_3D_computer_vision\Assignments\GEO1016_group10\Assignment1\A1_Calibration\A1_Calibration_Code\cmake-build-debug\CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object 3rd_party/glog/CMakeFiles/3rd_glog.dir/logging.cc.obj"
	cd C:\Users\noort\Documents\1\TU\year_4\GEO1016_3D_computer_vision\Assignments\GEO1016_group10\Assignment1\A1_Calibration\A1_Calibration_Code\cmake-build-debug\3rd_party\glog
	C:\PROGRA~2\MICROS~3\2019\BUILDT~1\VC\Tools\MSVC\1427~1.291\bin\Hostx86\x86\cl.exe @<<
 /nologo /TP $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) /FoCMakeFiles\3rd_glog.dir\logging.cc.obj /FdCMakeFiles\3rd_glog.dir\3rd_glog.pdb /FS -c C:\Users\noort\Documents\1\TU\year_4\GEO1016_3D_computer_vision\Assignments\GEO1016_group10\Assignment1\A1_Calibration\A1_Calibration_Code\3rd_party\glog\logging.cc
<<
	cd C:\Users\noort\Documents\1\TU\year_4\GEO1016_3D_computer_vision\Assignments\GEO1016_group10\Assignment1\A1_Calibration\A1_Calibration_Code\cmake-build-debug

3rd_party\glog\CMakeFiles\3rd_glog.dir\logging.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/3rd_glog.dir/logging.cc.i"
	cd C:\Users\noort\Documents\1\TU\year_4\GEO1016_3D_computer_vision\Assignments\GEO1016_group10\Assignment1\A1_Calibration\A1_Calibration_Code\cmake-build-debug\3rd_party\glog
	C:\PROGRA~2\MICROS~3\2019\BUILDT~1\VC\Tools\MSVC\1427~1.291\bin\Hostx86\x86\cl.exe > CMakeFiles\3rd_glog.dir\logging.cc.i @<<
 /nologo /TP $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E C:\Users\noort\Documents\1\TU\year_4\GEO1016_3D_computer_vision\Assignments\GEO1016_group10\Assignment1\A1_Calibration\A1_Calibration_Code\3rd_party\glog\logging.cc
<<
	cd C:\Users\noort\Documents\1\TU\year_4\GEO1016_3D_computer_vision\Assignments\GEO1016_group10\Assignment1\A1_Calibration\A1_Calibration_Code\cmake-build-debug

3rd_party\glog\CMakeFiles\3rd_glog.dir\logging.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/3rd_glog.dir/logging.cc.s"
	cd C:\Users\noort\Documents\1\TU\year_4\GEO1016_3D_computer_vision\Assignments\GEO1016_group10\Assignment1\A1_Calibration\A1_Calibration_Code\cmake-build-debug\3rd_party\glog
	C:\PROGRA~2\MICROS~3\2019\BUILDT~1\VC\Tools\MSVC\1427~1.291\bin\Hostx86\x86\cl.exe @<<
 /nologo /TP $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) /FoNUL /FAs /FaCMakeFiles\3rd_glog.dir\logging.cc.s /c C:\Users\noort\Documents\1\TU\year_4\GEO1016_3D_computer_vision\Assignments\GEO1016_group10\Assignment1\A1_Calibration\A1_Calibration_Code\3rd_party\glog\logging.cc
<<
	cd C:\Users\noort\Documents\1\TU\year_4\GEO1016_3D_computer_vision\Assignments\GEO1016_group10\Assignment1\A1_Calibration\A1_Calibration_Code\cmake-build-debug

# Object files for target 3rd_glog
3rd_glog_OBJECTS = \
"CMakeFiles\3rd_glog.dir\logging.cc.obj"

# External object files for target 3rd_glog
3rd_glog_EXTERNAL_OBJECTS =

lib\3rd_glog.lib: 3rd_party\glog\CMakeFiles\3rd_glog.dir\logging.cc.obj
lib\3rd_glog.lib: 3rd_party\glog\CMakeFiles\3rd_glog.dir\build.make
lib\3rd_glog.lib: 3rd_party\glog\CMakeFiles\3rd_glog.dir\objects1.rsp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=C:\Users\noort\Documents\1\TU\year_4\GEO1016_3D_computer_vision\Assignments\GEO1016_group10\Assignment1\A1_Calibration\A1_Calibration_Code\cmake-build-debug\CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX static library ..\..\lib\3rd_glog.lib"
	cd C:\Users\noort\Documents\1\TU\year_4\GEO1016_3D_computer_vision\Assignments\GEO1016_group10\Assignment1\A1_Calibration\A1_Calibration_Code\cmake-build-debug\3rd_party\glog
	$(CMAKE_COMMAND) -P CMakeFiles\3rd_glog.dir\cmake_clean_target.cmake
	cd C:\Users\noort\Documents\1\TU\year_4\GEO1016_3D_computer_vision\Assignments\GEO1016_group10\Assignment1\A1_Calibration\A1_Calibration_Code\cmake-build-debug
	cd C:\Users\noort\Documents\1\TU\year_4\GEO1016_3D_computer_vision\Assignments\GEO1016_group10\Assignment1\A1_Calibration\A1_Calibration_Code\cmake-build-debug\3rd_party\glog
	C:\PROGRA~2\MICROS~3\2019\BUILDT~1\VC\Tools\MSVC\1427~1.291\bin\Hostx86\x86\link.exe /lib /nologo /machine:X86 /out:..\..\lib\3rd_glog.lib @CMakeFiles\3rd_glog.dir\objects1.rsp 
	cd C:\Users\noort\Documents\1\TU\year_4\GEO1016_3D_computer_vision\Assignments\GEO1016_group10\Assignment1\A1_Calibration\A1_Calibration_Code\cmake-build-debug

# Rule to build all files generated by this target.
3rd_party\glog\CMakeFiles\3rd_glog.dir\build: lib\3rd_glog.lib

.PHONY : 3rd_party\glog\CMakeFiles\3rd_glog.dir\build

3rd_party\glog\CMakeFiles\3rd_glog.dir\clean:
	cd C:\Users\noort\Documents\1\TU\year_4\GEO1016_3D_computer_vision\Assignments\GEO1016_group10\Assignment1\A1_Calibration\A1_Calibration_Code\cmake-build-debug\3rd_party\glog
	$(CMAKE_COMMAND) -P CMakeFiles\3rd_glog.dir\cmake_clean.cmake
	cd C:\Users\noort\Documents\1\TU\year_4\GEO1016_3D_computer_vision\Assignments\GEO1016_group10\Assignment1\A1_Calibration\A1_Calibration_Code\cmake-build-debug
.PHONY : 3rd_party\glog\CMakeFiles\3rd_glog.dir\clean

3rd_party\glog\CMakeFiles\3rd_glog.dir\depend:
	$(CMAKE_COMMAND) -E cmake_depends "NMake Makefiles" C:\Users\noort\Documents\1\TU\year_4\GEO1016_3D_computer_vision\Assignments\GEO1016_group10\Assignment1\A1_Calibration\A1_Calibration_Code C:\Users\noort\Documents\1\TU\year_4\GEO1016_3D_computer_vision\Assignments\GEO1016_group10\Assignment1\A1_Calibration\A1_Calibration_Code\3rd_party\glog C:\Users\noort\Documents\1\TU\year_4\GEO1016_3D_computer_vision\Assignments\GEO1016_group10\Assignment1\A1_Calibration\A1_Calibration_Code\cmake-build-debug C:\Users\noort\Documents\1\TU\year_4\GEO1016_3D_computer_vision\Assignments\GEO1016_group10\Assignment1\A1_Calibration\A1_Calibration_Code\cmake-build-debug\3rd_party\glog C:\Users\noort\Documents\1\TU\year_4\GEO1016_3D_computer_vision\Assignments\GEO1016_group10\Assignment1\A1_Calibration\A1_Calibration_Code\cmake-build-debug\3rd_party\glog\CMakeFiles\3rd_glog.dir\DependInfo.cmake --color=$(COLOR)
.PHONY : 3rd_party\glog\CMakeFiles\3rd_glog.dir\depend

