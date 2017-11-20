#!/bin/bash
#PBS -A rck-371-aa
#PBS -l walltime=00:30:00
#PBS -l nodes=1:ppn=12
#PBS -q hb
#PBS -o outputfile
#PBS -e errorfile
#PBS -V
#PBS -N pz

# The executable to use
EXECUTABLE=@CMAKE_BINARY_DIR@/bin/test_integration_geometry

# Command line arguments
ARGV="extern_mesh/TEST_curved_2d_mixed"

# Specify the number of processor to run on (this should have correspondence with 'nodes' above)
N_PROCS="1"

# Specify the name of the file where output should be placed (leave empty to use stdout)
LOGFILE=""

# Specify the path to the mpi executable (mpiexec)
MPI_DIR=""

VALGRIND_OPTS="valgrind \
                 --track-origins=yes \
                 --leak-check=yes \
                 --num-callers=20 \
                 --suppressions=../external/valgrind/valgrind.supp \
              "
#                 --gen-suppressions=all \


# DO NOT MODIFY ANYTHING BELOW THIS LINE

#clear
${MPI_DIR}mpiexec -n ${N_PROCS} ${VALGRIND_OPTS} ${EXECUTABLE} ${ARGV} ${LOGFILE}
#${VALGRIND_OPTS} ${MPI_DIR}mpiexec -n ${N_PROCS} ${EXECUTABLE} ${ARGV} ${LOGFILE}
