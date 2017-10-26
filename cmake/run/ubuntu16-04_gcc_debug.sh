#!/bin/bash

TOP_DIR="${PWD}/../.."


# Modifiable parameters ****************************************************** #

# Must use mpich configured with --disable-checkpointing when running with valgrind otherwise there
# is a clash between the two programs. See [this][https://stackoverflow.com/a/37643501/5983549] SO
# answer.
export CMAKE_PREFIX_PATH=/home/pzwan/Applications/mpich/mpich-3.2/build

BUILD_DIR=${TOP_DIR}/build

CMAKE_BUILD_TYPE=Debug
TOOLCHAIN_FILE=gcc.cmake

# End Modifiable parameters ************************************************** #


mkdir -p ${BUILD_DIR} && cd ${BUILD_DIR}
cmake -D CMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE} \
      -D CMAKE_TOOLCHAIN_FILE=${TOP_DIR}/cmake/toolchains/${TOOLCHAIN_FILE} \
      ${TOP_DIR}
