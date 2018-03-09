#!/bin/bash

TOP_DIR="${PWD}"

cd ${TOP_DIR}/build_test_1D && make -j && make meshes
cd ${TOP_DIR}/build_test_2D && make -j && make meshes
cd ${TOP_DIR}/build_test_3D && make -j && make meshes
