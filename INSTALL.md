# Installation Instructions

The full list of software used in this project can be found on the
[code documentation page](https://codedocs.xyz/PhilipZwanenburg/DPGSolver/). The minimal list of programs necessary to
build/run the code is composed of the following:
- CMake
- MPICH
- Intel MKL
- GNU GSL
- PETSc
- Gmsh
- Python3

While it should be possible to install the majority of the packages using either `homebrew' (macOS) or `apt' (Ubuntu)
the following may need to be built manually:
- Intel MKL ([Download page](https://software.intel.com/en-us/mkl) or
[install using apt](https://software.intel.com/en-us/articles/installing-intel-free-libs-and-python-apt-repo))
- PETSc ([Download page](https://www.mcs.anl.gov/petsc/download/index.html))
- Gmsh ([Download page](http://gmsh.info/#Download))

The compressed files should be downloaded from the pages linked above and installed using the available gui or standard
configure, make, make install. Example detailed instructions for the PETSc installation are provided in the 
[PETSc README](external/petsc/) file.
