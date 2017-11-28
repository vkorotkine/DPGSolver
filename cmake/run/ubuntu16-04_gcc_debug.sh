#!/bin/bash

TOP_DIR="${PWD}/../.."


# Modifiable parameters ****************************************************** #

# Must use mpich configured with --disable-checkpointing when running with valgrind otherwise there
# is a clash between the two programs. See [this][https://stackoverflow.com/a/37643501/5983549] SO
# answer.
export CMAKE_PREFIX_PATH=/home/pzwan/Applications/mpich/mpich-3.2/build

export PETSC_DIR=/home/pzwan/Applications/petsc/petsc-3.8.0
export PETSC_ARCH=arch-linux-mkl-mpich-debug


BUILD_DIR=${TOP_DIR}/build

CMAKE_BUILD_TYPE=Debug
TOOLCHAIN_FILE=gcc.cmake

# End Modifiable parameters ************************************************** #

for dim in `seq 1 3`; do
	BUILD_DIR_D=${BUILD_DIR}_${dim}D
	mkdir -p ${BUILD_DIR_D} && cd ${BUILD_DIR_D}
	cmake -D CMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE} \
	      -D CMAKE_TOOLCHAIN_FILE=${TOP_DIR}/cmake/toolchains/${TOOLCHAIN_FILE} \
	      -D BUILD_DIM=${dim} \
	      ${TOP_DIR}
done
