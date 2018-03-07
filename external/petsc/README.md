# Installing PETSc

## Dependencies

The following packages must have previously been installed:
- [MPICH](https://www.mpich.org/);
- [Intel MKL](https://software.intel.com/en-us/mkl).
- Optionally: [ViennaCL](http://viennacl.sourceforge.net)

## Installation Instructions

- Clone PETSc in a suitable location using the command specified on the
  [PETSc download page][PETSc_download].
- From the PETSc ROOT directory, configure using the provided configure files. PETSC_ARCH will be
  set to the name of your configure file. Note that minor customization may be required for your
  system:

```sh
PETSC_ROOT$ ./relative_path_to_configure_file/PETSC_ARCH.py
```

- Follow the PETSc`make` instructions (which appear in a terminal message at the end of the previous
  stage):

```sh
PETSC_ROOT$ make BUILD_SPECIFIC_TEXT_HERE all
PETSC_ROOT$ make BUILD_SPECIFIC_TEXT_HERE test
```

The directory containing the libraries and include files for your architecture will now be present
in PETSC_ROOT.

## Options files

When required, the name of a PETSc options file from ```DPGROOT/external/petsc/options_files```
**without the `.txt` file extension** should be passed as argv[2].

[PETSc_download]: https://www.mcs.anl.gov/petsc/download/index.html
