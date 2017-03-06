# Additional Setup

Before running the code, be sure to execute:
$ make directories

Generate meshes using .geo files and place in appropriate folders created using make directories.

To run the code, cd to $(ROOT)/cases and execute using the script (.sh) files:
- quick.sh: Standard
- memcheck.sh: With valgrind enabled

Unit/Integration testing can be performed by uncommenting DTEST in the Makefile
