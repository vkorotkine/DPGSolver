#!/bin/bash
#PBS -A rck-371-aa
#PBS -l walltime=00:30:00
#PBS -l nodes=1:ppn=12
#PBS -q hb
#PBS -o outputfile
#PBS -e errorfile
#PBS -V
#PBS -N pz

# Executable and command line arguments
EXECUTABLE=@CMAKE_BINARY_DIR@/bin/test_integration_linearization
ARGV="navier_stokes/steady/taylor_couette/dg/TEST_NavierStokes_TaylorCouette_DG_TRI__ml0__p2"
ARGV="navier_stokes/steady/taylor_couette/dg/TEST_NavierStokes_TaylorCouette_DG_ParametricMixed2D__ml0__p1"
EXECUTABLE=@CMAKE_BINARY_DIR@/bin/test_integration_boundary_conditions
ARGV="boundary_conditions_navier_stokes_2d_noslip_general/adiabatic integration/boundary_conditions/TEST_NavierStokes_NSGA_2d__ml0__p2"

# Specify the number of processor to run on (this should have correspondence with 'nodes' above)
N_PROCS="1"

# Specify the name of the file where output should be placed (leave empty to use stdout)
LOGFILE=""

VALGRIND_OPTS="valgrind \
                 --track-origins=yes \
                 --leak-check=yes \
                 --num-callers=20 \
                 --suppressions=../external/valgrind/valgrind.supp \
              "
#                 --gen-suppressions=all \


# DO NOT MODIFY ANYTHING BELOW THIS LINE

@MPIEXEC@ -n ${N_PROCS} ${VALGRIND_OPTS} ${EXECUTABLE} ${ARGV} ${LOGFILE}
#${VALGRIND_OPTS} @MPIEXEC@ -n ${N_PROCS} ${EXECUTABLE} ${ARGV} ${LOGFILE}
