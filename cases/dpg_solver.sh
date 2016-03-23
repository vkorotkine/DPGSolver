#!/bin/bash
#PBS -A rck-371-aa
#PBS -l walltime=00:30:00
#PBS -l nodes=1:ppn=12
#PBS -q hb
#PBS -o outputfile
#PBS -e errorfile
#PBS -V
#PBS -N pz

LOGFILE="logfile"


TESTCASE="PeriodicVortex"
#TESTCASE="SupersonicVortex"


#USE_VALGRIND="0"
USE_VALGRIND="1"

case "${OSTYPE}" in
	*darwin*)
		CODE_DIR="/Users/philipzwanenburg/Desktop/Research_Codes"
		TOP_DIR="${CODE_DIR}/DPGC"
		MPI_DIR="${CODE_DIR}/Downloaded/petsc/petsc-3.6.3/arch-darwin-mpich-c-debug/bin/"
		N_PROCS="1"

		;;
	*linux*)
		TOP_DIR="/home/pzwan/Git/DPG"
		MPI_DIR=""
		#Make sure this is modified with "nodes" above
		N_PROCS="1"

		;;
esac

if [ "$USE_VALGRIND" = "1" ]; then
  VALGRIND_OPTS="valgrind --track-origins=yes --leak-check=yes"
else
  VALGRIND_OPTS=""
fi

clear
cd ${TOP_DIR}/cases

${MPI_DIR}mpiexec -n ${N_PROCS} ${VALGRIND_OPTS} ${TOP_DIR}/bin/DPGSolver.exe ${TESTCASE}
#${MPI_DIR}/mpiexec -n ${N_PROCS} ${VALGRIND_OPTS} ${TOP_DIR}/bin/DPGSolver.exe ${TESTCASE} > ${LOGFILE}


alias BEGINCOMMENT="if [ ]; then"
alias ENDCOMMENT="fi"

# Other options
BEGINCOMMENT

	valgrind --leak-check=yes
	valgrind --track-origins=yes --leak-check=yes

ENDCOMMENT
