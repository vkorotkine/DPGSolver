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
#include "Macros.h"
#include "S_DB.h"
#include "S_VOLUME.h"

#include "adaptation.h"
#include "update_VOLUMEs.h"
#include "implicit_VOLUME_info.h"
#include "implicit_FACET_info.h"
#include "finalize_LHS.h"
#include "output_to_paraview.h"

/*
 *	Purpose:
 *		Perform the implicit solve using Petsc's KSP object.
 *
 *	Comments:
 *		Using the residual as the initial guess for the iterative KSP solve resulted in divergence for the Poisson case.
 *		Chih-Hao mentioned that he never uses a non-zero initial guess and has not had problems. (ToBeModified)
 *
 *		Petsc's Cholesky solvers (direct and indirect) are much slower than the LU solvers (Petsc 3.6.3). ToBeModified
 *		
 *		Likely include a dynamic rtol value for KSPSetTolerances.
 *
 *	Notation:
 *
 *	References:
 */

void setup_KSP(Mat A, KSP ksp)
{
	// Initialize DB Parameters
	char *TestCase = DB.TestCase;

	// Standard datatypes
	char SolverType = 'i'; // Options: (i)terative, (d)irect

	// Petsc datatypes
	PC pc;

	KSPSetOperators(ksp,A,A);
	KSPSetTolerances(ksp,1e-10,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);
	KSPSetComputeSingularValues(ksp,PETSC_TRUE);

	KSPGetPC(ksp,&pc);
	if (strstr(TestCase,"Poisson")) {
		if (SolverType == 'i') {
			// Iterative Solve (Using Incomplete Cholesky)
			KSPSetInitialGuessNonzero(ksp,PETSC_TRUE);

			KSPSetType(ksp,KSPCG);

//			PCSetType(pc,PCICC);
			PCSetType(pc,PCILU);
			PCFactorSetLevels(pc,1); // Cannot use MatOrdering with 0 fill
			PCFactorSetMatOrderingType(pc,MATORDERINGRCM);
		} else {
			// Direct Solve (Using Cholesky Factorization)
			KSPSetType(ksp,KSPPREONLY);
//			PCSetType(pc,PCCHOLESKY);
			PCSetType(pc,PCLU);
		}
	} else {
		if (SolverType == 'i') {
			// Iterative Solve (Using ILU(1) with (R)everse (C)uthill-(M)cKee ordering)
			KSPSetInitialGuessNonzero(ksp,PETSC_TRUE);

			KSPSetType(ksp,KSPGMRES);
			KSPGMRESSetOrthogonalization(ksp,KSPGMRESModifiedGramSchmidtOrthogonalization);
//			KSPGMRESSetRestart(ksp,60); // Default: 30

			PCSetType(pc,PCILU);
			PCFactorSetLevels(pc,1); // Cannot use MatOrdering with 0 fill
			PCFactorSetMatOrderingType(pc,MATORDERINGRCM);
		} else {
			// Direct Solve (Using LU Factorization)
			KSPSetType(ksp,KSPPREONLY);
			PCSetType(pc,PCLU);
		}
	}
	KSPSetUp(ksp);
	PCSetUp(pc);
}

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

//	output_to_paraview("ZTest_Sol_Init");

	iteration = 0;
	maxRHS = 1.0; maxRHS0 = 1.0;
	while (maxRHS0/maxRHS < 1e10) {
		if (Adapt && iteration)
			mesh_update();

		// Build the RHS and LHS terms
		Mat                A = NULL;
		Vec                b = NULL, x = NULL;
		KSP                ksp;
//		PC                 pc;
		KSPConvergedReason reason;

		PetscInt *ix;

		printf("V");  implicit_VOLUME_info();
		printf("F");  implicit_FACET_info();
		printf("F "); maxRHS = finalize_LHS(&A,&b,&x,0);

		// Solve linear system
		printf("S");
		KSPCreate(MPI_COMM_WORLD,&ksp);
		setup_KSP(A,ksp);

		printf("S ");
		KSPSolve(ksp,b,x);
		KSPGetConvergedReason(ksp,&reason);
		KSPGetIterationNumber(ksp,&iteration_ksp);
//		KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD);

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

		printf("Iteration: %5d, KSP iterations (reason): %5d (%d), maxRHS (no MInv): % .3e\n",
		       iteration,iteration_ksp,reason,maxRHS);

		// Additional exit conditions
		if (maxRHS < 10*EPS && iteration) {
			printf("Exiting: maxRHS is below 10*EPS.\n");
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
