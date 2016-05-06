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


USE_VALGRIND="0"
#USE_VALGRIND="1"

KERNEL=$(uname -s)
case ${KERNEL} in
*darwin*)
	CODE_DIR="/Users/philipzwanenburg/Desktop/Research_Codes"
	TOP_DIR="${CODE_DIR}/DPGC"
	MPI_DIR="${CODE_DIR}/Downloaded/petsc/petsc-3.6.3/arch-darwin-mpich-c-debug/bin/"
	N_PROCS="1"

	;;
*Linux*)
	OS_RELEASE=$(uname -r)
	case ${OS_RELEASE} in
	*4.4.0-21-generic*)
		echo Linux
		echo ${OS_RELEASE}

		PROG_DIR="/home/philip/Desktop/research/programs"
		TOP_DIR="/home/philip/Desktop/research/codes/DPGSolver"
		MPI_DIR="${PROG_DIR}/petsc-3.7.0/arch-linux-c-/bin/"
		N_PROCS="1"

		;;
	*2.6.32-504.30.3.el6.x86_64*)
		PROG_DIR="/home/pzwan/programs"
		TOP_DIR="/home/pzwan/Git/DPG"
		MPI_DIR="${PROG_DIR}/petsc-3.6.3/arch-linux-mpich-c-opt/bin/"
		#Make sure this is modified with "nodes" above
		#N_PROCS="2"
		N_PROCS="1"

		;;
	esac

	;;
*linux*)
	echo linux not Linux? move guillimin stuff here

	;;
esac

case "${OSTYPE}" in
	*darwin*)
		CODE_DIR="/Users/philipzwanenburg/Desktop/Research_Codes"
		TOP_DIR="${CODE_DIR}/DPGC"
		MPI_DIR="${CODE_DIR}/Downloaded/petsc/petsc-3.6.3/arch-darwin-mpich-c-debug/bin/"
		N_PROCS="1"
		echo darwin

		;;
	*linux*)
    PROG_DIR="/home/pzwan/programs"
		TOP_DIR="/home/pzwan/Git/DPG"
		MPI_DIR="${PROG_DIR}/petsc-3.6.3/arch-linux-mpich-c-opt/bin/"
		#Make sure this is modified with "nodes" above
		#N_PROCS="2"
		N_PROCS="1"
		echo linux

		;;
esac
echo test

if [ "$USE_VALGRIND" = "1" ]; then
  VALGRIND_OPTS="valgrind --track-origins=yes --leak-check=yes"
  #VALGRIND_OPTS="valgrind --track-origins=yes --leak-check=yes --leak-check=full --show-reachable=yes"
else
  VALGRIND_OPTS=""
fi

#clear

echo test
echo ${OS}
echo "OS"

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
