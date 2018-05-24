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
ARGV="euler/TEST_Euler_PeriodicVortex_QUAD__ml0__p2"
ARGV="advection/gaussian_bump/dg/TEST_Advection_gaussian_bump2D petsc_options_gmres_default"

# Specify the number of processor to run on (this should have correspondence with 'nodes' above)
N_PROCS="1"

# Specify the name of the file where output should be placed (leave empty to use stdout)
LOGFILE=""

# Specify the path to the mpi executable (mpiexec)
MPI_DIR=""


# DO NOT MODIFY ANYTHING BELOW THIS LINE

#clear
EXECUTABLE=@CMAKE_BINARY_DIR@/bin/test_integration_convergence
ARGV="advection/gaussian_bump/dg/TEST_Advection_gaussian_bump1D petsc_options_gmres_default"
${MPI_DIR}mpiexec -n ${N_PROCS} ${EXECUTABLE} ${ARGV} ${LOGFILE}
EXECUTABLE=@CMAKE_BINARY_DIR@/bin/test_integration_convergence
ARGV="advection/gaussian_bump/dg/TEST_Advection_gaussian_bump2D petsc_options_gmres_default"
${MPI_DIR}mpiexec -n ${N_PROCS} ${EXECUTABLE} ${ARGV} ${LOGFILE}
EXECUTABLE=@CMAKE_BINARY_DIR@/bin/test_integration_convergence
ARGV="advection/gaussian_bump/dg/TEST_Advection_gaussian_bump3D petsc_options_gmres_default"
${MPI_DIR}mpiexec -n ${N_PROCS} ${EXECUTABLE} ${ARGV} ${LOGFILE}
