#!/usr/bin/env python2.7
configure_options = [
  '--with-mpi-dir=/usr/local/Cellar/mpich/3.2.1_2/',
  '--with-blas-lapack-dir=/Users/philip/Desktop/research_codes/mkl/mkl_2018_2/mkl',
  '--with-debugging=1',
  ]

if __name__ == '__main__':
  import sys,os
  sys.path.insert(0,os.path.abspath('config'))
  import configure
  configure.petsc_configure(configure_options)
