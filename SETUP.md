# Additional Setup

The configure/user_configure.mk file must be created before the code can be run. See the instructions in default.mk.

Before running the code, be sure execute the following in the $(ROOT) directory:
```sh
$ make directories
$ make meshes (see note below)
$ make
```

Note that making meshes currently requires that you add a 'user' name for yourself specifying the appropriate links to
the gmsh executable and the solver ROOT directory in 'meshfile_classes.py' as well as modifying the user in:
- MeshVariables_remove_duplicates.py
- MeshVariables_update.py
- generate_meshes.py

To run the code, cd to $(ROOT)/cases and execute using the script (.sh) files:
- quick.sh:    Standard
- memcheck.sh: With valgrind enabled

Parameters to be used for the run are specified in the '.ctrl' files located in $(ROOT)/cases/control_files.

Unit/Integration testing can be performed by uncommenting DTEST in the $(ROOT)/Makefile.
