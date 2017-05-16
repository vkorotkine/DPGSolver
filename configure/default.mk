# A customized copy of this file is required for each user. Please specify appropriate paths here and then save a copy
# of this file with the name 'user_configure.mk'

PROG_PATH  := 

PETSC_DIR  := 
PETSC_ARCH := 
METIS_DIR  := 
MKL_DIR    := 

CC := mpicc

# Parameters from Intal Math Kernel Library Link Line Advisor
# It was required that MKL be statically linked on OSX as the -Wl,--no-as-needed option is not supported by the OSX linker
MKL_LINKING := STATIC/DYNAMIC
MKL_INTERFACE_LAYER := 32/64
