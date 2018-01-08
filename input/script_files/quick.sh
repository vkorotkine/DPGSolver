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
ARGV="euler/supersonic_vortex/TEST_Euler_SupersonicVortex_DG_BlendedQUAD2D__ml0"
EXECUTABLE=@CMAKE_BINARY_DIR@/bin/test_integration_geometry
ARGV="extern_mesh/TEST_blended_2d_mixed"
EXECUTABLE=@CMAKE_BINARY_DIR@/bin/test_integration_convergence
ARGV="advection/default/dg/TEST_Advection_Default_DG_LINE__ml0 petsc_options_gmres_default"

# Specify the number of processor to run on (this should have correspondence with 'nodes' above)
N_PROCS="1"

# Specify the name of the file where output should be placed (leave empty to use stdout)
LOGFILE=""


# DO NOT MODIFY ANYTHING BELOW THIS LINE

@MPIEXEC@ -n ${N_PROCS} ${EXECUTABLE} ${ARGV} ${LOGFILE}
