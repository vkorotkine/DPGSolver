# Execute cmake using specific build parameters from '$ROOT/cmake/run' using, for example:
# $ sh gcc_debug.sh
#
# It is possible that linkage errors will occur when combining libraries in different ways as required files may be
# implicitly linked through other libraries because of CMake's treatment of "Cyclic Dependencies of Static Libraries".
#
# The minimum required version is required in order to properly find python3 without passing the cmake parameters:
# - PYTHON_INCLUDE_DIR, PYTHON_LIBRARY
#
# CMake Warnings which may be safely ignored when the conditions in braces () are satisfied:
# - "Manually-specified variables were not used by the project: CMAKE_TOOLCHAIN_FILE" (after the initial build)

cmake_minimum_required(VERSION 3.9)

message("Running cmake (version ${CMAKE_VERSION}) with build type: ${CMAKE_BUILD_TYPE}.\n")

project(DPGSolver)


set(DPGSolver_VERSION_MAJOR 1)
set(DPGSolver_VERSION_MINOR 0)
set(DPGSolver_VERSION "${DPGSolver_VERSION_MAJOR}.${DPGSolver_VERSION_MINOR}")

set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# Global includes
list	(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/modules"
	                          "${PROJECT_SOURCE_DIR}/cmake/modules/petsc_additional"
	)
include(functions)
include(ctest_testing)
include(multidimensional_builds)

# There was a problem on macOS (CMake 3.9.3) where enabling BUILD_SHARED_LIBS led to differing function addresses in
# different files.
option(BUILD_SHARED_LIBS "Build the shared library" OFF)

include_directories(src)

add_subdirectory(external)
add_subdirectory(src)
add_subdirectory(doc)
add_subdirectory(input)
