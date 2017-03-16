#!/bin/bash

GMSH="/home/philip/Desktop/research/programs/gmsh/gmsh-2.12.0-Linux/bin/gmsh"

${GMSH} TRI.geo -2 -setnumber MeshType 2
${GMSH} TRI.geo -2 -setnumber MeshType 3

#${GMSH} TRI.msh
