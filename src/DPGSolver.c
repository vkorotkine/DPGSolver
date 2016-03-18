#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include <petscksp.h>

#include "database.h"
#include "parameters.h"
#include "functions.h"

int main(int nargc, char **argv) {

  int MPIrank, MPIsize;

  // Start MPI and PETSC
  PetscInitialize(&nargc,&argv,PETSC_NULL,PETSC_NULL);
  MPI_Comm_size(MPI_COMM_WORLD,&MPIsize);
  MPI_Comm_rank(MPI_COMM_WORLD,&MPIrank);

  // Test memory leaks only from Petsc and MPI using valgrind
  //PetscFinalize(), exit(1);

  DB.MPIsize = MPIsize;
  DB.MPIrank = MPIrank;


  // Initialization
  Initialization(nargc,argv);

  // Preprocessing
  if (!DB.MPIrank) printf("Preprocessing:\n\n");

  if (!DB.MPIrank) printf("  Set up Parameters\n");
    SetupParameters();

  if (!DB.MPIrank) printf("  Set up Mesh\n");
    SetupMesh();
  



  // ONCE PETSC FUNCTIONS ARE USED, SWITCH TO PETSC-3.6

  MemoryFree();

  // End MPI and PETSC
  PetscFinalize();

  return 0;
}
