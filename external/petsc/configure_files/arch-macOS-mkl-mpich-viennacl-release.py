#!/usr/bin/env python2.7
configure_options = [
  '--with-mpi-dir=/Users/philip/Desktop/research_codes/mpich/mpich-3.2.1/mpich-install/',
  '--with-blas-lapack-dir=/Users/philip/Desktop/research_codes/mkl/mkl_2017_3/mkl',
  '--with-debugging=0',
  'COPTFLAGS=-O3 -march=native -mtune=native',
  'CXXOPTFLAGS=-O3 -march=native -mtune=native',
  'FOPTFLAGS=-O3 -march=native -mtune=native',
  '--with-viennacl=1',
  '--download-viennacl',
  '--with-openmp=1',
  'CFLAGS=-UVIENNACL_WITH_OPENCL -DVIENNACL_WITH_OPENMP',
  'CXXFLAGS=-UVIENNACL_WITH_OPENCL -DVIENNACL_WITH_OPENMP',
  ]

if __name__ == '__main__':
  import sys,os
  sys.path.insert(0,os.path.abspath('config'))
  import configure
  configure.petsc_configure(configure_options)
