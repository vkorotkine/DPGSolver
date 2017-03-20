# Additional Setup

Before running the code, be sure to execute:
$ make directories
$ make meshes

To run the code, cd to $(ROOT)/cases and execute using the script (.sh) files:
- quick.sh:    Standard
- memcheck.sh: With valgrind enabled

Parameters to be used for the run are specified in the '.ctrl' files situated in $(ROOT)/cases/control_files.

Unit/Integration testing can be performed by uncommenting DTEST in the ROOT Makefile.
