# Library and Program Installation Instructions

### Supported operating systems

Tested on OSX (10.12.5), Ubuntu (16.04)

### Detailed installation instructions

1. Install MPI and BLAS/LAPACK
- mpich2
- Intel MKL (tested with mkl_2017.3)

2. Install PETSc (tested with petsc-3.7.6)
- Create a custom configure file as in PETSC_ROOT/config/examples or as in the [example PETSc configure file](configure/example_PETSc_configure.py)
- Configure from PETSC_ROOT using the configure file as ./PATH_TO_PETSc_CONFIG/configure.py
- Follow make instructions for: make, make install, make test

3. Install ParMETIS (tested with parmetis-4.0.3)
- Download and extract .tar file
- In the root directory, modity the following parameters in the Makefile: cc = mpicc, cxx = mpicxx, debug = 0, BUILDDIR = build/name
- make config, make install
