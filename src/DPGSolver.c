#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include <petscksp.h>

#include "database.h"
#include "parameters.h"
#include "functions.h"

int main(int nargc, char **argv) {

  int MPIrank, MPIsize;
  struct S_DB DB;

  // Start MPI and PETSC
  PetscInitialize(&nargc,&argv,PETSC_NULL,PETSC_NULL);
  MPI_Comm_size(MPI_COMM_WORLD,&MPIsize);
  MPI_Comm_rank(MPI_COMM_WORLD,&MPIrank);

  DB.MPIsize = MPIsize;
  DB.MPIrank = MPIrank;

  // Initialization
  Initialization(nargc,argv);
  



  // ONCE PETSC FUNCTIONS ARE USED, SWITCH TO PETSC-3.6




  // End MPI and PETSC
  PetscFinalize();

  return 0;
}
