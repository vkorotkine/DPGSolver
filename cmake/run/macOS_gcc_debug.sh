#!/bin/bash

TOP_DIR="${PWD}/../.."


# Modifiable parameters ****************************************************** #

export MKLROOT=/Users/philip/Desktop/research_codes/mkl/mkl_2017_3/mkl

BUILD_DIR=${TOP_DIR}/build

CMAKE_BUILD_TYPE=Debug
TOOLCHAIN_FILE=gcc.cmake

# End Modifiable parameters ************************************************** #


mkdir -p ${BUILD_DIR} && cd ${BUILD_DIR}
cmake -D CMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE} \
      -D CMAKE_TOOLCHAIN_FILE=${TOP_DIR}/cmake/toolchains/${TOOLCHAIN_FILE} \
      ${TOP_DIR}
