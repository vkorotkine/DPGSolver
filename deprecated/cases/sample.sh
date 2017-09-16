#!/bin/bash
#PBS -A rck-371-aa
#PBS -l walltime=00:30:00
#PBS -l nodes=1:ppn=12
#PBS -q hb
#PBS -o outputfile
#PBS -e errorfile
#PBS -V
#PBS -N pz

# Specify the test case to run (Not used (but still required) when DTEST is enabled in the Makefile)
#TESTCASE="PeriodicVortex"
TESTCASE="SupersonicVortex"
#TESTCASE="dSphericalBump"
#TESTCASE="GaussianBump"
#TESTCASE="PolynomialBump"

# Specify the number of processor to run on (this should have correspondence with 'nodes' above)
N_PROCS="1"

# Specify the name of the file where output should be placed (leave empty to use stdout)
LOGFILE=""
#LOGFILE=" > logfile"

# Specify the path to the mpi executable (mpiexec)
MPI_DIR="/Users/philip/Desktop/research_codes/petsc/petsc-3.7.4/arch-osx-mpich-c-opt/bin/"
#MPI_DIR=""


# Specify whether the code should be run using valgrind and choose your options
USE_VALGRIND="0"
#USE_VALGRIND="1"

if [ "$USE_VALGRIND" = "1" ]; then
  VALGRIND_OPTS="valgrind --track-origins=yes --leak-check=yes"
else
  VALGRIND_OPTS=""
fi



# ADDITIONAL OPTIONS
alias BEGINCOMMENT="if [ ]; then"
alias ENDCOMMENT="fi"

# Various valgrind options
BEGINCOMMENT
	valgrind --track-origins=yes
	         --leak-check=yes
	         --leak-check=full
	         --show-reachable=yes
ENDCOMMENT



# DO NOT MODIFY ANYTHING BELOW THIS LINE

CURRENT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
EXEC_DIR="${CURRENT_DIR}/../bin"

clear

${MPI_DIR}mpiexec -n ${N_PROCS} ${VALGRIND_OPTS} ${EXEC_DIR}/DPGSolver.exe ${TESTCASE} ${LOGFILE}
