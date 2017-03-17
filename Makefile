# Makefile

# References
#   GNU Make manual
#   Miller(2008)-Recursive_Make_Considered_Harmful

# Notes
# - .RECIPEPREFIX == \t (TAB)

# Additional make targets:
# $ make directories
#
# $ make clean
# $ make clean_code
# $ make clean_test
# $ make clean_exec

# C standard
CSTD := -std=c99

# Options
OPTS := -O3
OPTS += -g -Wall -Wextra -Werror
OPTS += -DTEST


# Standard libraries (Math)
STD_LIB := -lm

# Machine dependent parameters
KERNEL  := $(shell uname -s)

LOCAL_INC := -I./include

# OSX
ifeq ($(KERNEL),Darwin)
#  PROG_PATH := /Users/philipzwanenburg/Desktop/Research_Codes/Downloaded
  PROG_PATH := /Users/philip/Desktop/research_codes

#  CC := $(PROG_PATH)/petsc/petsc-3.6.3/arch-darwin-mpich-c-debug/bin/mpicc -fopenmp -m64
#  CC := mpicc -fopenmp -m64

# There is a problem that the -fopenmp -m64 flag is not included in CC. Compiling fine for now => fix later. (ToBeDeleted)

  PETSC_DIR := $(PROG_PATH)/petsc/petsc-3.7.4
  PETSC_ARCH := arch-osx-mpich-c-opt
  CC := $(PETSC_DIR)/$(PETSC_ARCH)/bin/mpicc -fopenmp -m64

  # METIS_DIR := $(PROG_PATH)/parmetis/parmetis-4.0.3
  METIS_DIR := $(PROG_PATH)/parmetis/parmetis-4.0.3/build/opt

  METIS_INC      := -I$(METIS_DIR)/metis/include
  METIS_LDINC    := -L$(METIS_DIR)/libmetis -lmetis
  PARMETIS_INC   := -I$(METIS_DIR)/include
  PARMETIS_LDINC := -L$(METIS_DIR)/libparmetis -lparmetis

  MKL_DIR := $(PROG_PATH)/intel_MKL/mkl
  MKL_INC   := -I$(MKL_DIR)/include
  # MKL statically linked on OSX as the -Wl,--no-as-needed option is not supported by the OSX linker
  MKL_LDINC := $(MKL_DIR)/lib/libmkl_intel_lp64.a $(MKL_DIR)/lib/libmkl_core.a $(MKL_DIR)/lib/libmkl_sequential.a -lpthread
endif

# LINUX
ifeq ($(KERNEL),Linux)
  NODENAME := $(shell uname -n)
  ifeq ($(NODENAME),philip-Aspire-XC-605) #Home
    PROG_PATH := /home/philip/Desktop/research/programs

#    CC   := $(PROG_PATH)/petsc/petsc-3.7.0/arch-linux-c-/bin/mpicc -fopenmp -m64
    CC := mpicc -fopenmp -m64

#    PETSC_DIR := $(PROG_PATH)/petsc/petsc-3.7.0
#    PETSC_ARCH := arch-linux-c-
    PETSC_DIR := /home/philip/petsc/petsc-3.7.5
    PETSC_ARCH := linux-c-debug

    METIS_DIR      := $(PROG_PATH)/parmetis/parmetis-4.0.3/build/opt
    METIS_INC      := -I$(METIS_DIR)/metis/include
    METIS_LDINC    := -L$(METIS_DIR)/libmetis -lmetis
    PARMETIS_INC   := -I$(METIS_DIR)/include
    PARMETIS_LDINC := -L$(METIS_DIR)/libparmetis -lparmetis

    MKL_DIR := $(PROG_PATH)/intel/mkl
    MKL_INC := -I$(MKL_DIR)/include
#    MKL_LDINC := -Wl,--no-as-needed -L$(MKL_DIR)/lib/intel64 -lmkl_intel_lp64 -lmkl_core -lmkl_gnu_thread -ldl -lpthread -lgomp
    # Needed when petsc is not linked with MKL during installation
    MKL_LIBDIR := $(MKL_DIR)/lib/intel64
    MKL_LDINC := $(MKL_LIBDIR)/libmkl_intel_lp64.so $(MKL_LIBDIR)/libmkl_core.so $(MKL_LIBDIR)/libmkl_sequential.so -lpthread
  else #Guillimin
    PROG_PATH := /home/pzwan/programs

    CC   := $(PROG_PATH)/petsc/petsc-3.6.3/arch-linux-mpich-c-opt/bin/mpicc -fopenmp -m64

    PETSC_DIR := $(PROG_PATH)/petsc-3.6.3
    PETSC_ARCH := arch-linux-mpich-c-opt

    METIS_DIR := $(PROG_PATH)/parmetis-4.0.3/build/opt
    METIS_INC      := -I$(METIS_DIR)/metis/include
    METIS_LDINC    := -L$(METIS_DIR)/libmetis -lmetis
    PARMETIS_INC   := -I$(METIS_DIR)/include
    PARMETIS_LDINC := -L$(METIS_DIR)/libparmetis -lparmetis

    MKL_DIR := /software/compilers/Intel/2015-15.0/composer_xe_2015.0.090/mkl
    MKL_INC   := -I$(MKL_DIR)/include
    MKL_LDINC := -Wl,--no-as-needed -L$(MKL_DIR)/lib/intel64 -lmkl_intel_lp64 -lmkl_core -lmkl_gnu_thread -ldl -lpthread -lgomp
  endif
endif

# PETSC's 'variables' makefile
include $(PETSC_DIR)/lib/petsc/conf/variables
PETSC_INC := $(PETSC_CC_INCLUDES)

# Note: Parmetis must be linked before metis
#LIBS := $(STD_LIB) $(PETSC_LIB) $(PARMETIS_LDINC) $(METIS_LDINC) $(MKL_LDINC)
#INCS := $(LOCAL_INC) $(PETSC_INC) $(PARMETIS_INC) $(METIS_INC) $(MKL_INC)
LIBS := $(STD_LIB) $(MKL_LDINC) $(PETSC_LIB) $(PARMETIS_LDINC) $(METIS_LDINC)
INCS := $(LOCAL_INC) $(MKL_INC) $(PETSC_INC) $(PARMETIS_INC) $(METIS_INC)

EXECUTABLE := DPGSolver.exe

SRCDIR  := src
INCDIR  := include
OBJDIR  := obj
DEPDIR  := depend
EXECDIR := bin

nullstring :=
space := $(nullstring) # Single space
comma := ,

INC_DIRS := $(subst -I,,$(INCS))
INC_DIRS := $(subst $(space),:,$(INC_DIRS))
vpath %.c $(SRCDIR)
vpath %.h $(INC_DIRS)


EXECUTABLE := $(addprefix $(EXECDIR)/,$(EXECUTABLE))

SOURCES := $(wildcard $(SRCDIR)/*.c)
HEADERS := $(wildcard $(INCDIR)/*.h)
OBJECTS := $(SOURCES:$(SRCDIR)/%.c=$(OBJDIR)/%.o)
DEPENDS := $(SOURCES:$(SRCDIR)/%.c=$(DEPDIR)/%.d)

SHELL := /bin/bash
### Default goal + Additional required rules ###

# Compile executable file (Default goal)
$(EXECUTABLE) : $(OBJECTS)
	@echo
	@echo Creating/updating: $@
	@$(CC) -o $@ $(OPTS) $^ $(INCS) $(LIBS)
	@echo

# Include dependencies (Must be placed after default goal)
include $(DEPENDS)

# Create object and dependency files
$(OBJDIR)/%.o : %.c
	@echo Creating/updating: $@
	@$(CC) $(OPTS) $(CSTD) -c -o $@ $< $(INCS)

$(DEPDIR)/%.d : %.c
	@gcc -MM -MG $< > $@; # Use first prerequisite only as other prereqs are included from the existing %.d file
	@sed -i -e 's|.*:|$(OBJDIR)/$*.o $(DEPDIR)/$*.d:|' $@
	@echo Creating/updating: $@

# Additional dependencies needed for modified dynamic memory allocation
DYN_MEM_DEPS_NOPATH := S_ELEMENT.h S_VOLUME.h S_FACE.h
DYN_MEM_DEPS        := $(INCDIR)/S_ELEMENT.h $(INCDIR)/S_VOLUME.h $(INCDIR)/S_FACE.h

$(DYN_MEM_DEPS_NOPATH) : $(DYN_MEM_DEPS)
$(DYN_MEM_DEPS) : memory_constructors.c
	@echo
	@echo Updating memory_constructor dependencies.
	@echo
	@touch $(DYN_MEM_DEPS)


# Create directories if not present
$(OBJECTS): | $(OBJDIR)
$(OBJDIR):
	mkdir -p $(OBJDIR)

$(DEPENDS): | $(DEPDIR)
$(DEPDIR):
	mkdir -p $(DEPDIR)

$(EXECUTABLE): | $(EXECDIR)
$(EXECDIR):
	mkdir -p $(EXECDIR)


OUTPUT_LIST   := paraview errors results
TESTCASE_LIST := Poisson SupersonicVortex InviscidChannel SubsonicNozzle PrandtlMeyer
# To BeModified (Remove unneeded meshcases folders which may have been created previously)
MESHCASE_LIST := Ringleb \
                 dm1-Spherical_Section \
                 Ellipsoidal_Section \
                 Annular_Section \
                 Annular_Section/SupersonicVortex \
                 HoldenRamp \
                 GaussianBump \
                 NacaSymmetric \
                 EllipsoidalBump \
                 JoukowskiSymmetric
# ToBeModified (Remove unneeded meshtypes)
MESHTYPE_LIST := TRI CurvedTRI CurvedQUAD ToBeCurvedTRI ToBeCurvedQUAD Mixed CurvedMixed ToBeCurvedStructuredTRI

OUTPUT_LIST   := $(subst $(space),$(comma),$(OUTPUT_LIST))
TESTCASE_LIST := $(subst $(space),$(comma),$(TESTCASE_LIST))
MESHCASE_LIST := $(subst $(space),$(comma),$(MESHCASE_LIST))
MESHTYPE_LIST := $(subst $(space),$(comma),$(MESHTYPE_LIST))


### Additional Rules ###
.PHONY : directories
directories:
	mkdir -p cases/{$(OUTPUT_LIST)}
#	mkdir -p cases/{$(OUTPUT_LIST)}/{$(TESTCASE_LIST)}/{$(MESHTYPE_LIST)}
#	mkdir -p meshes/{$(MESHCASE_LIST)}
#	mkdir -p meshes/Test
	@echo

.PHONY : meshes

# Python compiler
PYTHONC := python3

MAIN_CONFIGURATIONS := Euler
TEST_CONFIGURATIONS := update_h linearization L2_proj Poisson Euler

MAIN_CONFIGURATIONS := $(addprefix main/,$(MAIN_CONFIGURATIONS))
TEST_CONFIGURATIONS := $(addprefix test/,$(TEST_CONFIGURATIONS))

CTRLDIR := cases/control_files

CONFIGURATIONS := $(addprefix $(CTRLDIR)/,$(MAIN_CONFIGURATIONS) $(TEST_CONFIGURATIONS))
CONTROL_FILES := $(shell find $(CTRLDIR) -name '*.ctrl')


meshes:
	@clear
	@clear
	@echo 
	@echo Creating MeshVariables file based on existing .ctrl files.
	@$(PYTHONC) python/update_MeshVariables.py $(CONFIGURATIONS)
	@echo 
#	$(MAKE) -C meshes test

$(CTRLDIR)/MeshVariables : $(CONTROL_FILES)




# Cleaning
.PHONY : clean clean_test clean_code clean_exec
clean:
	rm $(EXECUTABLE) $(OBJECTS) $(DEPENDS)

clean_test:
	rm $(OBJDIR)/test* $(DEPDIR)/test*

clean_code:
	find $(OBJDIR)/ -type f -not -name 'test*' -delete
	find $(DEPDIR)/ -type f -not -name 'test*' -delete

clean_exec:
	rm $(OBJDIR)/main.o $(DEPDIR)/main.d
