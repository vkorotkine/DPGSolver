#!/bin/bash

TOP_DIR="${PWD}"

cd ${TOP_DIR}/build_debug_1D && make -j && make meshes
cd ${TOP_DIR}/build_debug_2D && make -j && make meshes
cd ${TOP_DIR}/build_debug_3D && make -j && make meshes
