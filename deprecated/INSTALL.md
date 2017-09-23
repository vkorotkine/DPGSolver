# Detailed Installation Instructions

Last tested: 2017/05/18

### Supported Operating Systems

- macOS (tested on 10.12.5)
- Ubuntu (tested on 16.04)

### Installation Instructions

1. Install MPI and BLAS/LAPACK
	- MPICH
	- Intel MKL (tested with mkl_2017.3)
		- [Install using apt](https://software.intel.com/en-us/articles/installing-intel-free-libs-and-python-apt-repo).
		- Ensure that the 'Cluster support' package is also installed.
		- Consult the MKL [link line advisor](https://software.intel.com/en-us/articles/intel-mkl-link-line-advisor) for required information for [user_configure.mk](configure/user_configure.mk)

2. Install PETSc (tested with petsc-3.7.6)
	- Create a custom configure file as in PETSC_ROOT/config/examples or as in the [example PETSc configure file](configure/example_PETSc_configure.py). Note that PETSC_ARCH will be set to the name of your configure file.
	- Configure from PETSC_ROOT using the configure file as 
		```sh
		$ ./path_to_PETSc_configure_file/$(PETSC_ARCH).py
		```
	- Follow the PETSc make instructions:
		```sh
		$ make PETSC_DIR=/path_to_your_dir/ PETSC_ARCH=your_arch all
		$ make PETSC_DIR=/path_to_your_dir/ PETSC_ARCH=your_arch test
		```

3. Install ParMETIS (tested with parmetis-4.0.3)
	- Download and extract .tar file
	- In the root directory, modify the following parameters in the Makefile:
		```make
		cc = $(PATH_TO_MPICC)/mpicc
		cxx = $(PATH_TO_MPICXX)/mpicxx
		debug = 0 or 1 (choose one)
		prefix = full_path_to/BUILDDIR
		BUILDDIR = build/your_build_name
		```
	- Follow make instructions (elaborated in Install.txt/BUILD.txt):
		```sh
		$ make config
		$ make install
		```

4. Install gmsh
	- Add the path to the gmsh executable to $PATH (e.g. in .bash_profile add "PATH=$PATH:/Applications/Gmsh.app/Contents/MacOS/").
	- Check that this was successful by obtaining some output when typing
		```sh
		$ gmsh --version
		```

5. Install python3
	- Also install numpy.
