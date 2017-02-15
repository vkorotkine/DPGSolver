# Library and Program Installation Instructions

## Guillimin (Linux-x86_64)
**Using Petsc modules on guillimin resulted in errors during runtime, petsc was thus configured manually.**


### PETSc (petsc-3.6.3) - MPICH
0. Prerequisites (loaded modules):
  - cmake/2.8.12.2
  - MKL/11.2
1. Download and extract .tar file.
2. petsc configuration options:
  - A single build directory was made with debugging disabled.
  - -march=native had to be removed when compiling on OSX.
  - Linux: ./configure PETSC_ARCH=arch-linux-mpich-c-opt --download-mpich --with-blas-lapack-dir=/software/compilers/Intel/2015-15.0/composer_xe_2015.0.090/mkl --with-debugging=0 COPTFLAGS='-O3 -march=native -mtune=native' CXXOPTFLAGS='-O3 -march=native -mtune=native' FOPTFLAGS='-O3 -march=native -mtune=native'
3. Follow make instructions for installation and testing.

### ParMETIS (parmetis-4.0.3)
0. Prerequisites:
  - MPI (either OpenMPI or mpich (used below), link appropriately)
1. Download and extract .tar file
2. In the top directory, modify the following parameters in the Makefile:
  - debug  = 1
  - prefix = /home/pzwan/programs/parmetis-4.0.3/build/opt/
  - cc     = /home/pzwan/programs/petsc-3.6.3/arch-linux-mpich-c-opt/bin/mpicc
  - cxx    = /home/pzwan/programs/petsc-3.6.3/arch-linux-mpich-c-opt/bin/mpicxx
  - BUILDDIR = build/opt
3. Execute make commands
  - make config
  - make install


## Personal Computer (OSX 10.10)
**Using MPICH resulted in fewer valgrind memory leaks.**


### PETSC (petsc-3.6.3) - MPICH
0. Prerequisites:
  - Intel MKL
1. Download and extract .tar file.
2. PETSC Configuration options:
  - It is suggested to make two build directories, with debugging enabled/disabled.
  - I was unable to build PETSC using COPTFLAGS='-O3 -march=native -mtune=native' on OSX 10.10
  - ./configure PETSC_ARCH=arch-osx-mpich-c-opt --download-mpich --with-blas-lapack-dir=/Users/philip/Desktop/research_codes/intel_MKL/compilers_and_libraries_2016.4.210/mac/mkl --with-debugging=0 COPTFLAGS='-O3 -mtune=native' CXXOPTFLAGS='-O3 -mtune=native' FOPTFLAGS='-O3 -mtune=native'
3. Follow make instructions for installation and testing.

### ParMETIS (parmetis-4.0.3)
0. Prerequisites:
  - CMake
  - MPI (either OpenMPI or mpich (used below), link appropriately)
1. Download and extract .tar file.
2. In the top directory, modify the following parameters in the Makefile:
  - debug  = 0
  - prefix = /Users/philip/Desktop/research_codes/parmetis/parmetis-4.0.3/build/opt/
  - cc     = /Users/philip/Desktop/research_codes/petsc/petsc-3.7.4/arch-osx-mpich-c-opt/bin/mpicc
  - cxx    = /Users/philip/Desktop/research_codes/petsc/petsc-3.7.4/arch-osx-mpich-c-opt/bin/mpicxx
  - BUILDDIR = build/opt
3. Execute make commands.
  - make config
  - make install
