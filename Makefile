# Makefile
# See the GNU Make manual for guidelines.

# Instructions:
# 	Place "Makefile" in the top directory
# 	Place source code "*.c" and "*.h" files in "src" folder
#
# Additional Options:
# 	make clean


# C standard and compiler
CSTD = -std=c99

# Options
#OPTS = -O3
OPTS = -g -Wall -Wextra -O3
#OPTS = -g -Wall -Wextra -O3 -DTEST

# Standard libraries (Math)
STD_LIB = -lm

# Machine dependent parameters
KERNEL  := $(shell uname -s)

#all:
#	@echo $(KERNEL)

# OSX
ifeq ($(KERNEL),Darwin)
  PROG_PATH = /Users/philipzwanenburg/Desktop/Research_Codes/Downloaded

  CC := ${PROG_PATH}/petsc/petsc-3.6.3/arch-darwin-mpich-c-debug/bin/mpicc -fopenmp -m64
#  CC = ${PROG_PATH}/petsc/petsc-3.6.3/arch-darwin-mpich-c-opt/bin/mpicc -fopenmp -m64
#  CC = mpicc -fopenmp -m64

# There is a problem that the -fopenmp -m64 flag is not included in CC. Compiling fine for now => fix later. (ToBeDeleted)

#  PETSC_DIR = ${PROG_PATH}/petsc/petsc-3.2-p7 (ToBeDeleted)
  PETSC_DIR = ${PROG_PATH}/petsc/petsc-3.6.3
  PETSC_ARCH = arch-darwin-mpich-c-debug
#  PETSC_ARCH = arch-darwin-mpich-c-opt
#  PETSC_ARCH = arch-darwin-c-opt

  # METIS_DIR = ${PROG_PATH}/parmetis/parmetis-4.0.3
  METIS_DIR = ${PROG_PATH}/parmetis_mpich/parmetis-4.0.3/build/debug
#  METIS_DIR = ${PROG_PATH}/parmetis_mpich/parmetis-4.0.3/build/opt

  METIS_INC      = -I${METIS_DIR}/metis/include
  METIS_LDINC    = -L${METIS_DIR}/libmetis -lmetis
  PARMETIS_INC   = -I${METIS_DIR}/include
  PARMETIS_LDINC = -L${METIS_DIR}/libparmetis -lparmetis

  MKL_DIR   = ${PROG_PATH}/intel/mkl
  MKL_INC = -I${MKL_DIR}/include
  # MKL statically linked on OSX as the -Wl,--no-as-needed option is not supported by the OSX linker
  MKL_LDINC = ${MKL_DIR}/lib/libmkl_intel_lp64.a ${MKL_DIR}/lib/libmkl_core.a ${MKL_DIR}/lib/libmkl_sequential.a -lpthread
endif

# LINUX
ifeq ($(KERNEL),Linux)
  OS_RELEASE := $(shell uname -r)
  ifeq ($(OS_RELEASE),4.4.0-21-generic) #Home
    PROG_PATH = /home/philip/Desktop/research/programs

    CC   = ${PROG_PATH}/petsc/petsc-3.7.0/arch-linux-c-/bin/mpicc -fopenmp -m64

    PETSC_DIR = ${PROG_PATH}/petsc/petsc-3.7.0
    PETSC_ARCH = arch-linux-c-

    METIS_DIR = ${PROG_PATH}/parmetis/parmetis-4.0.3/build/opt
    METIS_INC      = -I${METIS_DIR}/metis/include
    METIS_LDINC    = -L${METIS_DIR}/libmetis -lmetis
    PARMETIS_INC   = -I${METIS_DIR}/include
    PARMETIS_LDINC = -L${METIS_DIR}/libparmetis -lparmetis

    MKL_DIR = ${PROG_PATH}/intel/mkl
    MKL_INC = -I${MKL_DIR}/include
    MKL_LDINC = -Wl,--no-as-needed -L$(MKL_DIR)/lib/intel64 -lmkl_intel_lp64 -lmkl_core -lmkl_gnu_thread -ldl -lpthread -lgomp
  else #Guillimin
    PROG_PATH = /home/pzwan/programs

    CC   = ${PROG_PATH}/petsc/petsc-3.6.3/arch-linux-mpich-c-opt/bin/mpicc -fopenmp -m64

    PETSC_DIR = ${PROG_PATH}/petsc-3.6.3
    PETSC_ARCH = arch-linux-mpich-c-opt

    METIS_DIR = ${PROG_PATH}/parmetis-4.0.3/build/opt
    METIS_INC      = -I${METIS_DIR}/metis/include
    METIS_LDINC    = -L${METIS_DIR}/libmetis -lmetis
    PARMETIS_INC   = -I${METIS_DIR}/include
    PARMETIS_LDINC = -L${METIS_DIR}/libparmetis -lparmetis

    MKL_DIR = /software/compilers/Intel/2015-15.0/composer_xe_2015.0.090/mkl
    MKL_INC   = -I$(MKL_DIR)/include
    MKL_LDINC = -Wl,--no-as-needed -L$(MKL_DIR)/lib/intel64 -lmkl_intel_lp64 -lmkl_core -lmkl_gnu_thread -ldl -lpthread -lgomp
  endif
endif

PETSC_INC = -I PETSC_DIR/include ${PETSC_CC_INCLUDES}
# Run PETSC's 'variables' makefile
include ${PETSC_DIR}/lib/petsc/conf/variables

# missing LIBPATH, DEFINES (as compared to Brian's makefile)
# Parmetis must be linked before metis
#LIBS = $(STD_LIB) $(PETSC_INC) $(PETSC_LIB) $(PARMETIS_INC) $(PARMETIS_LDINC) $(METIS_INC) $(METIS_LDINC) $(MKL_INC) $(MKL_LDINC)
LIBS = $(STD_LIB) $(PETSC_LIB) $(PARMETIS_LDINC) $(METIS_LDINC) $(MKL_LDINC)
INCS = $(PETSC_INC) $(PARMETIS_INC) $(METIS_INC) $(MKL_INC)

EXECUTABLE = DPGSolver.exe

SRCDIR  = src
OBJDIR  = obj
EXECDIR = bin

# OBJECTS    := $(addprefix $(OBJDIR)/,$(OBJECTS))
EXECUTABLE := $(addprefix $(EXECDIR)/,$(EXECUTABLE))

SOURCES = $(wildcard $(SRCDIR)/*.c)
HEADERS = $(wildcard $(SRCDIR)/*.h)
#OBJECTS = $(patsubst $(SRCDIR)/%.c,$(OBJDIR)/%.o,$(wildcard $(SRCDIR)/*.c))
OBJECTS = $(SOURCES:$(SRCDIR)/%.c=$(OBJDIR)/%.o)


# Formatting for "rules" in makefiles is as follows:
# 	target ... : prerequisites ...
# 		recipe
# 		...
#
# 	Note that a "tab" character must be placed before each line of the recipe
#
# Automatic variables:
#
# 	$@: target
# 	$^: prerequisites
# 	$<: first prerequisite only (desired if excluding headers for example)
#
# Compile commands:
# 	-c:       no linking is done
# 	-o <arg>: output to <arg>

### Default goal + Additional required rules ###

# Compile executable file (Default goal)
$(EXECUTABLE) : $(OBJECTS)
	@echo 
	@echo Building executable file
	@echo
	$(CC) -o $@ $(OPTS) $^ $(INCS) $(LIBS)

# Create objects
# Still need to figure out how to include header dependencies.
# See Miller(2008)-Recursive_Make_Considered_Harmful
$(OBJECTS) : $(OBJDIR)/%.o : $(SRCDIR)/%.c
	@echo 
	@echo Building/updating object files
	@echo
	$(CC) $(OPTS) $(CSTD) -c -o $@ $< $(INCS)

# Create directories if not present
$(OBJECTS): | $(OBJDIR)
$(OBJDIR):
	mkdir $(OBJDIR)

$(EXECUTABLE): | $(EXECDIR)
$(EXECDIR):
	mkdir $(EXECDIR)


### Additional Rules ###

# Make clean ($ make clean)
# ".PHONY" is included here to specify that clean is a phony target

.PHONY : clean
clean:
	rm $(EXECUTABLE) $(OBJECTS)
