#!/bin/bash

TOP_DIR="${PWD}"

cd ${TOP_DIR}/build_1D && make && make meshes
cd ${TOP_DIR}/build_2D && make && make meshes
cd ${TOP_DIR}/build_3D && make && make meshes
