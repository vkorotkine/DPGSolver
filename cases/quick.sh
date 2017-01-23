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

#TESTCASE="PeriodicVortex"
TESTCASE="SupersonicVortex"
#TESTCASE="dSphericalBump"
#TESTCASE="GaussianBump"
#TESTCASE="PolynomialBump"


USE_VALGRIND="0"
#USE_VALGRIND="1"

KERNEL=$(uname -s)
case ${KERNEL} in
*arwin*)
	CODE_DIR="/Users/philip/Desktop/"
	PROG_DIR="/Users/philip/Desktop/research_codes"
	TOP_DIR="${CODE_DIR}/DPGSolver"
	MPI_DIR="${PROG_DIR}/petsc/petsc-3.7.4/arch-osx-mpich-c-opt/bin/"
	N_PROCS="1"

	;;
*inux*)
	NODENAME=$(uname -n)
	case ${NODENAME} in
	*philip-Aspire-XC-605*)
#		echo Linux
#		echo ${OS_RELEASE}

		PROG_DIR="/home/philip/Desktop/research/programs"
		TOP_DIR="/home/philip/Desktop/research/codes/DPGSolver"
#		MPI_DIR="${PROG_DIR}/petsc/petsc-3.7.0/arch-linux-c-/bin/"
		MPI_DIR=
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
esac

if [ "$USE_VALGRIND" = "1" ]; then
  VALGRIND_OPTS="valgrind --track-origins=yes --leak-check=yes"
  #VALGRIND_OPTS="valgrind --track-origins=yes --leak-check=yes --leak-check=full --show-reachable=yes"
else
  VALGRIND_OPTS=""
fi

clear

cd ${TOP_DIR}/cases

${MPI_DIR}mpiexec -n ${N_PROCS} ${VALGRIND_OPTS} ${TOP_DIR}/bin/DPGSolver.exe ${TESTCASE}
#${MPI_DIR}mpiexec -n ${N_PROCS} -ksp_monitor_true_residual ${VALGRIND_OPTS} ${TOP_DIR}/bin/DPGSolver.exe ${TESTCASE}
#${MPI_DIR}/mpiexec -n ${N_PROCS} ${VALGRIND_OPTS} ${TOP_DIR}/bin/DPGSolver.exe ${TESTCASE} > ${LOGFILE}


alias BEGINCOMMENT="if [ ]; then"
alias ENDCOMMENT="fi"

# Other options
BEGINCOMMENT

	valgrind --leak-check=yes
	valgrind --track-origins=yes --leak-check=yes

ENDCOMMENT
