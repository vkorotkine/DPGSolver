#!/bin/bash

TOP_DIR="${PWD}/../.."


# Modifiable parameters ****************************************************** #

export MKLROOT=/opt/intel/mkl
export PETSC_DIR=/Users/manmeetbhabra/Research_Codes/petsc-3.8.4
export PETSC_ARCH=arch-macOS-mkl-mpich-debug

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
