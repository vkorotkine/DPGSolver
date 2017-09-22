#!/usr/bin/env python
configure_options = [
  '--with-cc=mpicc',
  '--with-cxx=mpicxx',
  '--with-fc=mpif90',
  '--with-blas-lapack-dir=/opt/intel/mkl',
  '--with-debugging=1',
  ]

if __name__ == '__main__':
  import sys,os
  sys.path.insert(0,os.path.abspath('config'))
  import configure
  configure.petsc_configure(configure_options)
