#!/usr/bin/env python
configure_options = [
  '--with-mpi-dir=/home/pzwan/Applications/mpich/mpich-3.2/build/',
  '--with-blas-lapack-dir=/opt/intel/mkl',
  '--with-debugging=0',
  'COPTFLAGS=-O3 -march=native -mtune=native',
  'CXXOPTFLAGS=-O3 -march=native -mtune=native',
  'FOPTFLAGS=-O3 -march=native -mtune=native',
  ]

if __name__ == '__main__':
  import sys,os
  sys.path.insert(0,os.path.abspath('config'))
  import configure
  configure.petsc_configure(configure_options)
