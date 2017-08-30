# Execute cmake using specific build parameters from '$ROOT/cmake/run' using, for example:
# $ sh gcc_debug.sh
#
# CMake Warnings which may be safely ignored when the conditions in braces () are satisfied:
# - "Manually-specified variables were not used by the project: CMAKE_TOOLCHAIN_FILE" (after the initial build)

cmake_minimum_required(VERSION 3.5.1)

message("Running cmake with build type: ${CMAKE_BUILD_TYPE}.\n")

project(DPGSolver)


set(DPGSolver_VERSION_MAJOR 1)
set(DPGSolver_VERSION_MINOR 0)
set(DPGSolver_VERSION "${DPGSolver_VERSION_MAJOR}.${DPGSolver_VERSION_MINOR}")


# Global includes
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/modules")
include(functions)

# Include subdirectory CMakeLists.txt files
include_directories(src)

add_subdirectory(src)
add_subdirectory(doc)
add_subdirectory(meshes)
