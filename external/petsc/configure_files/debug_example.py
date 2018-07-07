#!/usr/bin/env python

# Make a copy of this file and rename it to the desired name for PETSC_ARCH.
configure_options = [
  '--with-mpi-dir=/usr/',
  '--with-blas-lapack-dir=/opt/intel/mkl',
  '--with-debugging=1',
  ]

if __name__ == '__main__':
  import sys,os
  sys.path.insert(0,os.path.abspath('config'))
  import configure
  configure.petsc_configure(configure_options)
