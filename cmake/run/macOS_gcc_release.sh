#!/bin/bash

DIR="${BASH_SOURCE%/*}"
if [[ ! -d "$DIR" ]]; then DIR="$PWD"; fi

PARAM_FILE="$DIR/parameters_release.sh"
PARAM_EXAMPLE="$DIR/parameters_example.sh"
if [ ! -f  ${PARAM_FILE} ]; then
	echo "Parameter file (${PARAM_FILE}) not found."
	echo "Please create it based on the example file (${PARAM_EXAMPLE})."
else
	. ${PARAM_FILE}

	TOP_DIR="${PWD}/../.."
	BUILD_DIR=${TOP_DIR}/build

	CMAKE_BUILD_TYPE=Release
	TOOLCHAIN_FILE=gcc.cmake

	for dim in `seq 1 3`; do
		BUILD_DIR_D=${BUILD_DIR}_${dim}D
		mkdir -p ${BUILD_DIR_D} && cd ${BUILD_DIR_D}
		cmake -D CMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE} \
			-D CMAKE_TOOLCHAIN_FILE=${TOP_DIR}/cmake/toolchains/${TOOLCHAIN_FILE} \
			-D BUILD_DIM=${dim} \
			${TOP_DIR}
	done
fi
