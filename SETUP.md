# Additional Set Up

An out of source build must be created (by CMake) using the provided bash scripts:
```sh
ROOT/cmake/run$ sh gcc_debug.sh
```

which creates a build director of ROOT/build by default. The following `make` commands are then available in the build
directory:
	- make meshes (Generate meshes)
	- make doc    (Generate documentation)

The code can then be run from using the sample bash scripts in BUILD/script_files.


# Additional Set Up (Deprecated Version)

## Configuring the Code

The 'configure/user_configure.mk' file must be created before the code can be run. See the instructions in [default.mk](configure/default.mk).

Before running the code, be sure to execute the following in the $(ROOT) directory:
```sh
$ make directories
$ make meshes
$ make
```

## Running the Code

1. cd to $ROOT/cases,
2. copy the [sample.sh](cases/sample.sh) script to 'your_script_file_name.sh' and modify the the parameters according to your build
3. run using
```sh
$ sh your_script_file_name.sh
```

Parameters to be used for the run are specified in the '.ctrl' files located in $(ROOT)/cases/control_files.

Unit/Integration testing can be performed by uncommenting DTEST in the [ROOT Makefile](Makefile).

## Generating Documentation

Copy the [default Doxygen configure file](configure/default_doxygen.cfg) into the [doc](doc/) directory after modifying
the 'STRIP_FROM_PATH' parameter. By default, documentation is not generated for static functions.

Generate the documentation by running Doxygen:
```{sh}
$ doxygen name_of_your_config.cfg
```
