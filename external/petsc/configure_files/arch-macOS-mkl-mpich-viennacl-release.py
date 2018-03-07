#!/usr/bin/env python2.7
configure_options = [
  '--with-mpi-dir=/usr/local/Cellar/mpich/3.2_3/',
  '--with-blas-lapack-dir=/Users/philip/Desktop/research_codes/mkl/mkl_2017_3/mkl',
  '--with-debugging=0',
  'COPTFLAGS=-O3 -march=native -mtune=native',
  'CXXOPTFLAGS=-O3 -march=native -mtune=native',
  'FOPTFLAGS=-O3 -march=native -mtune=native',
  '--with-viennacl=1',
  '--with-viennacl-dir=/Users/philip/Desktop/research_codes/viennacl/ViennaCL-1.7.1/build',
#  '--download-viennacl',
#  '--with-opencl-include=/System/Library/Frameworks/OpenCL.framework/Versions/Current/Headers',
#  '--with-opencl-lib=-framework OpenCL',
  ]

if __name__ == '__main__':
  import sys,os
  sys.path.insert(0,os.path.abspath('config'))
  import configure
  configure.petsc_configure(configure_options)
