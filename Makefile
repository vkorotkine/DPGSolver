# Makefile for most code functionality

# References
#   GNU Make manual
#   Miller(2008)-Recursive_Make_Considered_Harmful

# Notes
# - .RECIPEPREFIX == \t (TAB)

# Additional make targets:
# directories
# meshes
#
# clean
# clean_code
# clean_test
# clean_exec
# clean_empty (removes empty directories)

# C standard
CSTD := -std=c11

# Options
OPTS := -O3
OPTS += -g -Wall -Wextra -Werror
OPTS += -DTEST


# Standard libraries (Math)
STD_LIB := -lm


DPG_ROOT := $(shell pwd)
ifeq (,$(wildcard $(DPG_ROOT)/configure/user_configure.mk))
   $(error The 'configure/user_configure.mk' file is not present. Please create it)
endif
include configure/user_configure.mk


LOCAL_INC := -I./include

ifeq ($(MKL_INTERFACE_LAYER),32)
    MKL_INC := 
else ifeq ($(MKL_INTERFACE_LAYER),64)
    MKL_INC := -DMKL_ILP64 -m64
endif
MKL_INC += -I$(MKL_DIR)/include


include $(PETSC_DIR)/lib/petsc/conf/variables
PETSC_INC := $(PETSC_CC_INCLUDES)

METIS_INC      := -I$(METIS_DIR)/metis/include
METIS_LDINC    := -L$(METIS_DIR)/libmetis -lmetis
PARMETIS_INC   := -I$(METIS_DIR)/include
PARMETIS_LDINC := -L$(METIS_DIR)/libparmetis -lparmetis

# Note: Parmetis must be linked before metis
LIBS := $(STD_LIB) $(MKL_LDINC) $(PETSC_LIB) $(PARMETIS_LDINC) $(METIS_LDINC)
INCS := $(LOCAL_INC) $(MKL_INC) $(PETSC_INC) $(PARMETIS_INC) $(METIS_INC)


EXECUTABLE := DPGSolver.exe

SRCDIR  := src
INCDIR  := include
OBJDIR  := obj
DEPDIR  := depend
EXECDIR := bin
CTRLDIR := cases/control_files
MESHDIR := meshes

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
TESTCASE_LIST := Advection_Default \
                 Advection_Peterson \
                 Poisson \
                 Euler_PeriodicVortex \
                 Euler_PeriodicVortex_Stationary \
                 Euler_SupersonicVortex \
                 Euler_InternalSubsonic_GaussianBump \
                 NavierStokes_TaylorCouette \
                 NavierStokes_PlaneCouette

MESHTYPE_LIST := LINE \
                 TRI CurvedTRI ToBeCurvedTRI \
                 QUAD CurvedQUAD ToBeCurvedQUAD \
                 TET ToBeCurvedTET \
                 HEX ToBeCurvedHEX \
                 WEDGE ToBeCurvedWEDGE \
                 PYR ToBeCurvedPYR \
                 MIXED2D CurvedMIXED2D ToBeCurvedMIXED2D \

OUTPUT_LIST   := $(subst $(space),$(comma),$(OUTPUT_LIST))
TESTCASE_LIST := $(subst $(space),$(comma),$(TESTCASE_LIST))
MESHTYPE_LIST := $(subst $(space),$(comma),$(MESHTYPE_LIST))


### Additional Rules ###
.PHONY : directories
directories:
	@echo
	@echo Creating directories if not present
	mkdir -p cases/{$(OUTPUT_LIST)}/{$(TESTCASE_LIST)}/{$(MESHTYPE_LIST)}
	@echo


# Python compiler
PYTHONC := python3

#MAIN_CONFIGURATIONS := Euler
MAIN_CONFIGURATIONS := $(nullstring)
TEST_CONFIGURATIONS := update_h L2_proj_p L2_proj_h Advection Poisson Euler NavierStokes

MAIN_CONFIGURATIONS := $(addprefix main/,$(MAIN_CONFIGURATIONS))
TEST_CONFIGURATIONS := $(addprefix test/,$(TEST_CONFIGURATIONS))

CONFIGURATIONS := $(addprefix $(CTRLDIR)/,$(MAIN_CONFIGURATIONS) $(TEST_CONFIGURATIONS))
CONTROL_FILES  := $(shell find $(CTRLDIR) -name '*.ctrl')
MESHVARIABLES  := $(MESHDIR)/MeshVariables $(MESHDIR)/MeshVariables_python


.PHONY : meshes
meshes:
	@echo
	$(MAKE) mesh_vars_and_deps
	$(MAKE) -C meshes meshes_all
	@echo


.PHONY : mesh_vars_and_deps
mesh_vars_and_deps : $(MESHVARIABLES)
$(MESHVARIABLES) : $(CONTROL_FILES)
	@echo
	@echo Creating MeshVariables file based on existing .ctrl files.
	cd python && $(PYTHONC) MeshVariables_update.py $(CONFIGURATIONS)
	@cd python && $(PYTHONC) MeshVariables_remove_duplicates.py
	cd python/peterson_mesh/ && $(PYTHONC) MeshVariables_update.py
	@echo



# Cleaning
.PHONY : clean clean_test clean_code clean_exec clean_empty
clean:
	rm $(EXECUTABLE) $(OBJECTS) $(DEPENDS)

clean_test:
	rm $(OBJDIR)/test* $(DEPDIR)/test*

clean_code:
	find $(OBJDIR)/ -type f -not -name 'test*' -delete
	find $(DEPDIR)/ -type f -not -name 'test*' -delete

clean_exec:
	rm $(OBJDIR)/main.o $(DEPDIR)/main.d

clean_empty:
	find . -type d -empty -delete
