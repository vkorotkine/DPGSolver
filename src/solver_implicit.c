// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

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
#include "S_ELEMENT.h"
#include "S_VOLUME.h"
#include "Test.h"

#include "adaptation.h"
#include "update_VOLUMEs.h"
#include "implicit_GradW.h"
#include "implicit_VOLUME_info.h"
#include "implicit_FACE_info.h"
#include "finalize_LHS.h"
#include "output_to_paraview.h"
#include "element_functions.h"
#include "matrix_functions.h"
#include "variable_functions.h"
#include "solver_explicit.h"

#include "array_print.h"

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
 *		Adytia recommended trying PCSetType(pc,PCGAMG) (Does not required SetLevels or SetMatOrderingType). May need to
 *		specify that blocks are arranged in sizes of Nvar*Neq*NvnS*NvnS. (ToBeModified)
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
	KSPSetTolerances(ksp,1e-15,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);
//	KSPSetTolerances(ksp,1e-10,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);
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
			KSPGMRESSetRestart(ksp,60); // Default: 30

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

static void compute_underRelax(struct S_VOLUME *VOLUME, const double *dWhat, double *alpha)
{
	/*
	 *	Purpose:
	 *		Compute underRelaxation scaling based on maintaining physically realistic density and pressure and limiting
	 *		their change by a certain percentage.
	 *
	 *	Comments:
	 *		As this only checks values at certain nodes within the element, this may still result in unphysical updates.
	 *		In case of this occurence, perhaps try to use a basis with a convex hull property so that the max and min in
	 *		the element can be known.
	 */

	double       maxChange = 0.8;

	// Initialize DB Parameters
	unsigned int d    = DB.d,
	             Nvar = DB.Nvar;

	// Standard datatypes
	unsigned int P, NvnS;
	double       *What, *ChiS_vS, *W0, *dW, *U0, *U, *dU;

	struct S_ELEMENT *ELEMENT;

	P    = VOLUME->P;
	NvnS = VOLUME->NvnS;
	What = VOLUME->What;

	ELEMENT = get_ELEMENT_type(VOLUME->type);

	ChiS_vS = ELEMENT->ChiS_vS[P][P][0];

	W0 = malloc(NvnS*Nvar * sizeof *W0); // free
	dW = malloc(NvnS*Nvar * sizeof *dW); // free

	U0 = malloc(NvnS*Nvar * sizeof *U0); // free
	U  = malloc(NvnS*Nvar * sizeof *U);  // free
	dU = malloc(NvnS*Nvar * sizeof *dU); // free

	mm_d(CBCM,CBT,CBNT,NvnS,Nvar,NvnS,1.0,0.0,ChiS_vS, What,W0);
	mm_d(CBCM,CBT,CBNT,NvnS,Nvar,NvnS,1.0,0.0,ChiS_vS,dWhat,dW);

	convert_variables(W0,U0,d,d,NvnS,1,'c','p');
	convert_variables(dW,dU,d,d,NvnS,1,'c','p');

	double *rho0, *drho, *p0, *dp, alphaO;

	rho0 = &U0[0]; p0 = &U0[(Nvar-1)*NvnS];
	drho = &dU[0]; dp = &dU[(Nvar-1)*NvnS];

	unsigned int flag[2] = {1,1};
	alphaO = 2.0*(*alpha);
//	alphaO = 1.0;
	while (flag[0] || flag[1]) {
		unsigned int n;

		alphaO /= 2.0;

		// rho
		flag[0] = 0;

		for (n = 0; n < NvnS; n++) {
			double rho = rho0[n]+alphaO*drho[n];
			if (rho <= 0.0 || rho/rho0[n] > (1+maxChange) || rho/rho0[n] < (1-maxChange)) {
//				printf("%d % .3e % .3e % .3e % .3e\n",n,alphaO,rho,rho0[n],drho[n]);
				flag[0] = 1;
				break;
			}
		}

		// p
		flag[1] = 0;

		for (n = 0; n < NvnS; n++) {
			if (flag[0])
				break;

			double p = p0[n]+alphaO*dp[n];
//printf("Out: %d %d % .3e % .3e % .3e\n",n,flag[0],p,p0[n],dp[n]);
			if (p <= 0.0 || p/p0[n] > (1+maxChange) || p/p0[n] < (1-maxChange)) {
				flag[1] = 1;
				break;
			}
		}

		if (alphaO < EPS) {
			printf("%d %d\n",flag[0],flag[1]);
			printf("Potential problem: Under relaxation driven to 0.\n");//, EXIT_MSG;
			break;
		}

	}
	if (alphaO < EPS)
		alphaO = 1e5*EPS;

	*alpha = alphaO;

	free(W0);
	free(dW);
	free(U0);
	free(U);
	free(dU);
}

void solver_implicit_linear_system(Mat *A, Vec *b, Vec *x, KSP *ksp, unsigned int const iteration,
                                   bool const PrintEnabled)
{
	if (PrintEnabled) { printf("F "); }
	double const maxRHS = finalize_LHS(A,b,x,0);

	bool const Output_A_MATLAB = 0;
	if (Output_A_MATLAB) {
		// Used for outputting matrix to file in matlab sparse format
		PetscViewer viewer;
		PetscViewerASCIIOpen(PETSC_COMM_WORLD,"mat.output", &viewer);
		PetscViewerPushFormat(viewer,PETSC_VIEWER_ASCII_MATLAB);
		MatView(*A,viewer);
		PetscViewerDestroy(&viewer);
	}

	// Setup Preconditioner
	if (PrintEnabled) { printf("S"); }
	KSPCreate(MPI_COMM_WORLD,ksp);
	setup_KSP(*A,*ksp);

	// Solve the system of equations
	if (PrintEnabled) { printf("S "); }
	KSPSolve(*ksp,*b,*x);

	KSPConvergedReason reason;
	KSPGetConvergedReason(*ksp,&reason);

	PetscInt iteration_ksp;
	KSPGetIterationNumber(*ksp,&iteration_ksp);

	PetscReal emax, emin;
	KSPComputeExtremeSingularValues(*ksp,&emax,&emin);

	// Display solver progress
	if (PrintEnabled) {
		printf("Iteration: %5d, KSP iterations (cond, reason): %5d (% .3e, %d), maxRHS (no MInv): % .3e\n",
		       iteration,iteration_ksp,emax/emin,reason,maxRHS);
	}
}

void solver_implicit_update_What(Vec x) {
	unsigned int Nvar = DB.Nvar;

	for (struct S_VOLUME *VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next) {
		unsigned int const IndA = VOLUME->IndA,
		                   NvnS = VOLUME->NvnS,
		                   Ndof = NvnS*Nvar;

		double *const What  = VOLUME->What,
		       *const dWhat = malloc(Ndof * sizeof *dWhat); // free

		PetscInt *const ix = malloc(Ndof * sizeof *ix); // free
		for (size_t i = 0; i < Ndof; i++)
			ix[i] = IndA+i;

		VecGetValues(x,Ndof,ix,dWhat);
		free(ix);

		for (size_t i = 0; i < Ndof; i++)
			What[i] += dWhat[i];
		free(dWhat);
	}
}

void solver_implicit(bool const PrintEnabled)
{
	// Initialize DB Parameters
	unsigned int OutputInterval = DB.OutputInterval,
	             Nvar           = DB.Nvar,
	             Adapt          = DB.Adapt;

	unsigned int PrintTesting = 1;

	// Standard datatypes
	char         *dummyPtr_c[2], *string, *fNameOut;
	unsigned int i, iMax, iteration, IndA, NvnS;
	int          iteration_ksp;
	double       maxRHS0, maxRHS, *What, *dWhat;

	struct S_VOLUME *VOLUME;

	for (i = 0; i < 2; i++)
		dummyPtr_c[i] = malloc(STRLEN_MIN * sizeof *dummyPtr_c[i]); // free

	fNameOut = malloc(STRLEN_MAX * sizeof *fNameOut); // free
	string   = malloc(STRLEN_MIN * sizeof *string);   // free

	update_VOLUME_Ops();
	update_VOLUME_finalize();

	output_to_paraview("ZTest_Sol_Init");

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

		if (PrintEnabled && DB.Viscous) { printf("G"); } implicit_GradW();

		if (PrintEnabled) { printf("V");  } implicit_VOLUME_info();
		if (PrintEnabled) { printf("F");  } implicit_FACE_info();
		if (PrintEnabled) { printf("F "); } maxRHS = finalize_LHS(&A,&b,&x,0);

		// Used for outputting matrix to file in matlab sparse format
//		PetscViewer viewer;
//		PetscViewerASCIIOpen(PETSC_COMM_WORLD,"mat.output", &viewer);
//		PetscViewerPushFormat(viewer,PETSC_VIEWER_ASCII_MATLAB);
//		MatView(A,viewer);
//		PetscViewerDestroy(&viewer);

		// Solve linear system
		if (PrintEnabled) { printf("S"); }
		KSPCreate(MPI_COMM_WORLD,&ksp);
		setup_KSP(A,ksp);

		if (PrintEnabled) { printf("S "); }
		KSPSolve(ksp,b,x);
		KSPGetConvergedReason(ksp,&reason);
		KSPGetIterationNumber(ksp,&iteration_ksp);

		PetscReal emax, emin;
		KSPComputeExtremeSingularValues(ksp,&emax,&emin);
//		KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD);

		// Display solver progress
		if (!iteration)
			maxRHS0 = maxRHS;

		if (PrintEnabled) {
			printf("Iteration: %5d, KSP iterations (cond, reason): %5d (% .3e, %d), maxRHS (no MInv): % .3e\n",
			       iteration,iteration_ksp,emax/emin,reason,maxRHS);
		}

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

			double alpha = 1.0;
if (iteration < 3)
			compute_underRelax(VOLUME,dWhat,&alpha);
			//printf("% .3e\n",alpha);

			for (i = 0; i < iMax; i++) {
				if (fabs(dWhat[i]) <= EPS)
					What++;
				else
					(*What++) += alpha*dWhat[i];
			}
			free(dWhat);

			enforce_positivity_highorder(VOLUME);
		}

		KSPDestroy(&ksp);
		finalize_ksp(&A,&b,&x,2);

		// Output to paraview
		if (PrintTesting && (iteration % OutputInterval == 0 || iteration < 3)) {
			sprintf(dummyPtr_c[1],"%d",iteration);
			strcpy(dummyPtr_c[0],"SolStart");
			strcat(dummyPtr_c[0],dummyPtr_c[1]);
			output_to_paraview(dummyPtr_c[0]);
		}

		// Additional exit conditions
		if (maxRHS < 10*EPS && iteration) {
			if (PrintEnabled) { printf("Exiting: maxRHS is below 10*EPS.\n"); }
			break;
		}

		// hp adaptation
//		if (Adapt)
		if (0&&Adapt)
			adapt_hp();

		iteration++;
	}

	// Output to paraview
	if (TestDB.ML <= 1 || (TestDB.PGlobal == 1) || (TestDB.PGlobal+TestDB.ML) <= 8) {
		strcpy(fNameOut,"SolFinal_");
		sprintf(string,"%dD_",DB.d);   strcat(fNameOut,string);
		                               strcat(fNameOut,DB.MeshType);
		if (DB.Adapt == ADAPT_0) {
			sprintf(string,"_ML%d",DB.ML); strcat(fNameOut,string);
			sprintf(string,"P%d_",DB.PGlobal); strcat(fNameOut,string);
		} else {
			sprintf(string,"_ML%d",TestDB.ML); strcat(fNameOut,string);
			sprintf(string,"P%d_",TestDB.PGlobal); strcat(fNameOut,string);
		}
		output_to_paraview(fNameOut);
	}

	for (i = 0; i < 2; i++)
		free(dummyPtr_c[i]);

	free(fNameOut);
	free(string);
}
