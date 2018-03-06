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
EXECUTABLE=@CMAKE_BINARY_DIR@/bin/test_integration_convergence
ARGV="RESULTS_navier_stokes_joukowski_full_dg_2D petsc_options_gmres_tol_1e-2"
ARGV="RESULTS_euler_supersonic_vortex_dg_2d petsc_options_gmres_tol_1e-2"

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
