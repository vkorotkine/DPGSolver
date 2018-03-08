#!/usr/bin/env python

# CURRENTLY NOT WORKING! PROBLEM WITH DOUBLE PRECISION.
configure_options = [
  '--with-mpi-dir=/home/pzwan/Applications/mpich/mpich-3.2/build/',
  '--with-blas-lapack-dir=/opt/intel/mkl',
  '--with-debugging=0',
  'COPTFLAGS=-O3 -march=native -mtune=native',
  'CXXOPTFLAGS=-O3 -march=native -mtune=native',
  'FOPTFLAGS=-O3 -march=native -mtune=native',
  '--download-viennacl=yes',
#  '--with-viennacl=1',
#  '--with-viennacl-dir=/home/pzwan/Applications/viennacl/ViennaCL-1.7.1/build',
  '--with-opencl=1',
  ]

if __name__ == '__main__':
  import sys,os
  sys.path.insert(0,os.path.abspath('config'))
  import configure
  configure.petsc_configure(configure_options)
