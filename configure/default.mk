# A customized copy of this file is required for each user. Please specify appropriate paths here and then save a copy
# of this file with the name 'user_configure.mk'

PROG_PATH  := 

PETSC_DIR  := 
PETSC_ARCH := 
METIS_DIR  := 
MKL_DIR    := 

CC := mpicc -m64

# Parameters from Intal Math Kernel Library Link Line Advisor
MKL_LINKING := DYNAMIC
MKL_INTERFACE_LAYER := 32 # Potentially update to 64 in future

LIB_DIR   := $(MKL_DIR)/lib
MKL_LDINC :=  -L$(LIB_DIR)/intel64 -Wl,-no-as-needed -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl
