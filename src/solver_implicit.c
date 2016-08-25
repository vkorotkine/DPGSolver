// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

#include "solver_implicit.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "petscvec.h"
#include "petscmat.h"
#include "petscksp.h"

#include "Parameters.h"
#include "S_DB.h"
#include "S_VOLUME.h"

#include "adaptation.h"
#include "update_VOLUMEs.h"
#include "implicit_VOLUME_info.h"
#include "implicit_FACET_info.h"
#include "finalize_LHS.h"
#include "output_to_paraview.h"

#include "Macros.h"
#include "array_print.h"

/*
 *	Purpose:
 *		Perform the implicit solve using Petsc's KSP object.
 *
 *	Comments:
 *
 *	Notation:
 *
 *	References:
 */

void solver_implicit(void)
{
	// Initialize DB Parameters
	unsigned int OutputInterval = DB.OutputInterval,
	             Nvar           = DB.Nvar,
	             Adapt          = DB.Adapt;

	unsigned int PrintTesting = 1;

	// Standard datatypes
	char         *dummyPtr_c[2];
	unsigned int i, iMax, iteration, IndA, NvnS;
	int          iteration_ksp;
	double       maxRHS0, maxRHS, *What, *dWhat;

	struct S_VOLUME *VOLUME;

	for (i = 0; i < 2; i++)
		dummyPtr_c[i] = malloc(STRLEN_MIN * sizeof *dummyPtr_c[i]); // free

	update_VOLUME_finalize();

output_to_paraview("ZTest_Sol_Init");

	iteration = 0;
	maxRHS = 1.0; maxRHS0 = 1.0;
	while (maxRHS0/maxRHS < 1e10) {
		if (Adapt && iteration)
			mesh_update();

		// Build the RHS and LHS terms
		Mat A = NULL;
		Vec b = NULL, x = NULL;
		KSP ksp;
//		PC  pc;

		PetscInt *ix;

		printf("V");  implicit_VOLUME_info();
		printf("F");  implicit_FACET_info();
		printf("F "); maxRHS = finalize_LHS(&A,&b,&x,0);

		// Solve linear system
		KSPCreate(MPI_COMM_WORLD,&ksp);

		KSPSetOperators(ksp,A,A);//,DIFFERENT_NONZERO_PATTERN); ToBeDeleted (This was included in Brian's code)
// Potentially modify preconditioner matrix (ToBeDeleted)
//		KSPSetInitialGuessNonzero(ksp,PETSC_TRUE);
//		KSPSetTolerances(ksp,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);
		KSPSetTolerances(ksp,1e-15,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);

//KSPGetPC(ksp,&pc);
//PCSetType(pc,PCJACOBI);
//PCSetType(pc,PCILU);

//		KSPSetTolerances(ksp,0.0,0.0,1e250,100);


// Modify parameters here (ToBeDeleted)

//		KSPSetFromOptions(ksp);
		KSPSetUp(ksp);

		KSPSolve(ksp,b,x);
//seems to be diverging...

PetscBool flg  = PETSC_FALSE;
PetscOptionsGetBool(NULL,"-ksp_reason",&flg,NULL);
printf("flg: %d\n",flg);
KSPConvergedReason reason;
KSPGetConvergedReason(ksp,&reason);
PetscPrintf(PETSC_COMM_WORLD,"KSPConvergedReason: %D\n", reason);

		KSPGetIterationNumber(ksp,&iteration_ksp);
printf("%d\n",iteration_ksp);

		// Update What
		for (VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next) {
			IndA = VOLUME->IndA;
			NvnS = VOLUME->NvnS;
			What = VOLUME->What;

			iMax = NvnS*Nvar;
			ix    = malloc(iMax * sizeof *ix);    // free
			dWhat = malloc(iMax * sizeof *dWhat); // free
			for (i = 0; i < iMax; i++)
				ix[i] = IndA+i;

			VecGetValues(x,iMax,ix,dWhat);
			free(ix);
if (1||iteration == 100) {
printf("\n");
array_print_d(NvnS,Nvar,What,'C');
array_print_d(NvnS,Nvar,dWhat,'C');
EXIT_MSG;
}

			for (i = 0; i < iMax; i++)
				(*What++) += dWhat[i];
			free(dWhat);
		}

		KSPDestroy(&ksp);
		finalize_ksp(&A,&b,&x,2);

		// Output to paraview
		if (PrintTesting && (iteration % OutputInterval == 0 || iteration < 5)) {
			sprintf(dummyPtr_c[1],"%d",iteration);
			strcpy(dummyPtr_c[0],"SolStart");
			strcat(dummyPtr_c[0],dummyPtr_c[1]);
			output_to_paraview(dummyPtr_c[0]);
		}

		// Display solver progress
		if (!iteration)
			maxRHS0 = maxRHS;

		printf("Iteration: %8d, KSP iterations: %8d, maxRHS (no MInv): % .3e\n",iteration,iteration_ksp,maxRHS);

		// Additional exit conditions
		if (maxRHS < 8e-14 && iteration > 2) {
			printf("Exiting: maxRHS is below 8e-14.\n");
			break;
		}

		// hp adaptation
		if (Adapt)
			adapt_hp();

		iteration++;
	}

	for (i = 0; i < 2; i++)
		free(dummyPtr_c[i]);
}
