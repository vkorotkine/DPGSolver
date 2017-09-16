# A customized copy of this file is required for each user. Please specify appropriate paths here and then save a copy
# of this file with the name 'user_configure.mk'

PROG_PATH  := /home/pzwan/Applications

PETSC_DIR  := $(PROG_PATH)/petsc
PETSC_ARCH := arch-linux-mkl-debug
METIS_DIR  := $(PROG_PATH)/parmetis/parmetis-4.0.3/build/arch-linux-debug
MKL_DIR    := /opt/intel/mkl

CC := mpicc -m64

# Parameters from Intel Math Kernel Library Link Line Advisor
MKL_LINKING := DYNAMIC
MKL_INTERFACE_LAYER := 32 # Potentially update to 64 in future

LIB_DIR   := $(MKL_DIR)/lib
MKL_LDINC :=  -L$(LIB_DIR)/intel64 -Wl,-no-as-needed -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl
