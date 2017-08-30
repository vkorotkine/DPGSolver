#!/bin/bash

TOP_DIR="${PWD}/../.."


# Modifiable parameters ****************************************************** #

BUILD_DIR=${TOP_DIR}/build

CMAKE_BUILD_TYPE=Debug
TOOLCHAIN_FILE=gcc.cmake

# End Modifiable parameters ************************************************** #


mkdir -p ${BUILD_DIR} && cd ${BUILD_DIR}
cmake -D CMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE} \
      -D PYTHON_INCLUDE_DIR=/usr/include/python3.5 \
      -D PYTHON_LIBRARY=/usr/lib/python3.5/config-3.5m-x86_64-linux-gnu/libpython3.5.so \
      -D CMAKE_TOOLCHAIN_FILE=${TOP_DIR}/cmake/toolchains/${TOOLCHAIN_FILE} \
      ${TOP_DIR}
