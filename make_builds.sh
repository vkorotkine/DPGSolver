#!/bin/bash

TOP_DIR="${PWD}"

cd ${TOP_DIR}/build_1D && make -j && make meshes
cd ${TOP_DIR}/build_2D && make -j && make meshes
cd ${TOP_DIR}/build_3D && make -j && make meshes
