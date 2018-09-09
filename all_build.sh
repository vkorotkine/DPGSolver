#!/bin/bash

TOP_DIR="${PWD}"

cd ${TOP_DIR}/cmake/run && ./debug.sh
cd ${TOP_DIR}/cmake/run && ./release.sh

