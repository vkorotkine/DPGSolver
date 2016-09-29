// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

#include "test_integration_linearization.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <complex.h>

#include "petscmat.h"

#include "Parameters.h"
#include "Macros.h"
#include "S_DB.h"
#include "S_VOLUME.h"
#include "S_FACET.h"
#include "Test.h"

#include "test_code_integration.h"
#include "test_support.h"
#include "array_norm.h"
#include "array_print.h"
#include "array_free.h"
#include "solver_poisson_c.h"
#include "explicit_VOLUME_info_c.h"
#include "implicit_VOLUME_info.h"
#include "explicit_FACET_info_c.h"
#include "implicit_FACET_info.h"
#include "finalize_RHS_c.h"
#include "finalize_LHS.h"

/*
 *	Purpose:
 *		Test correctness of implementation of code linearization.
 *
 *	Comments:
 *		By default, the complete linearization is checked here. However, linearizations of the individual contributions
 *		listed below may be checked separately (See the commented code in 'test_integration_linearization'):
 *			0) complete check (default)
 *			1) diagonal VOLUME contributions
 *			2) diagonal FACET contributions
 *			3) off-diagonal FACET contributions
 *		Further, the complete linearization can be checked either through the assembly of the individual contributions
 *		or more simply by using the assembled RHS directly.
 *
 *	Notation:
 *
 *	References:
 */

void compute_A_cs(Mat *A, Vec *b, Vec *x, const unsigned int assemble_type)
{
	if (!assemble_type) {
		compute_A_cs(A,b,x,1);
		compute_A_cs(A,b,x,2);
		compute_A_cs(A,b,x,3);

		finalize_Mat(A,1);
		return;
	}

	// Initialize DB Paramters
	char         *TestCase = DB.TestCase;
	unsigned int Nvar      = DB.Nvar;

	// Standard datatypes
	unsigned int   i, j, iMax, jMax, side, NvnS[2], IndA[2], nnz_d;
	double         h;
	double complex *RHS_c;

	struct S_VOLUME *VOLUME, *VOLUME2;
	struct S_FACET  *FACET;

	PetscInt    *m, *n;
	PetscScalar *vv;

	h = EPS*EPS;

	if (*A == NULL)
		initialize_KSP(A,b,x);

	for (VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next) {
		NvnS[0] = VOLUME->NvnS;

		// Note: Initialize with zeros for linear cases.
		if (strstr(TestCase,"Poisson")) {
			if (VOLUME->uhat_c)
				free(VOLUME->uhat_c);

			VOLUME->uhat_c = calloc(NvnS[0]*Nvar , sizeof *(VOLUME->uhat_c));
		} else {
			if (VOLUME->What_c)
				free(VOLUME->What_c);

			VOLUME->What_c = malloc(NvnS[0]*Nvar * sizeof *(VOLUME->What_c));

			for (i = 0, iMax = NvnS[0]*Nvar; i < iMax; i++)
				VOLUME->What_c[i] = VOLUME->What[i];
		}
	}

	for (VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next) {
		IndA[0] = VOLUME->IndA;
		NvnS[0] = VOLUME->NvnS;
		for (i = 0, iMax = NvnS[0]*Nvar; i < iMax; i++) {
			if (strstr(TestCase,"Poisson")) {
				VOLUME->uhat_c[i] += h*I;

				compute_qhat_VOLUME_c();
				compute_qhat_FACET_c();
				finalize_qhat_c();
				switch (assemble_type) {
				default: // 0
					compute_uhat_VOLUME_c();
					compute_uhat_FACET_c();
					break;
				case 1:
					compute_uhat_VOLUME_c();
					break;
				case 2:
				case 3:
					compute_uhat_FACET_c();
					break;
				}
			} else {
				VOLUME->What_c[i] += h*I;

				switch (assemble_type) {
				default: // 0
					explicit_VOLUME_info_c();
					explicit_FACET_info_c();
					break;
				case 1:
					explicit_VOLUME_info_c();
					break;
				case 2:
				case 3:
					explicit_FACET_info_c();
					break;
				}
			}

			if (assemble_type == 1) {
				for (VOLUME2 = DB.VOLUME; VOLUME2; VOLUME2 = VOLUME2->next) {
					if (VOLUME->indexg != VOLUME2->indexg)
						continue;

					IndA[1] = VOLUME2->IndA;
					NvnS[1] = VOLUME2->NvnS;
					nnz_d   = VOLUME2->nnz_d;

					m  = malloc(NvnS[1]*Nvar * sizeof *m);  // free
					n  = malloc(NvnS[1]*Nvar * sizeof *n);  // free
					vv = malloc(NvnS[1]*Nvar * sizeof *vv); // free

					RHS_c = VOLUME2->RHS_c;
					for (j = 0, jMax = NvnS[1]*Nvar; j < jMax; j++) {
						m[j]  = IndA[1]+j;
						n[j]  = IndA[0]+i;
						vv[j] = cimag(RHS_c[j])/h;
					}

					MatSetValues(*A,nnz_d,m,1,n,vv,ADD_VALUES);
					free(m); free(n); free(vv);
				}
			}

			if (assemble_type == 2 || assemble_type == 3) {
				for (FACET = DB.FACET; FACET; FACET = FACET->next) {
				for (side = 0; side < 2; side++) {
					if (assemble_type == 2) {
						if (side == 0) {
							VOLUME2 = FACET->VIn;
							RHS_c = FACET->RHSIn_c;
						} else {
							if (FACET->Boundary)
								continue;
							VOLUME2 = FACET->VOut;
							RHS_c = FACET->RHSOut_c;
						}
						if (VOLUME->indexg != VOLUME2->indexg)
							continue;

						IndA[1] = VOLUME2->IndA;
						NvnS[1] = VOLUME2->NvnS;
						nnz_d   = VOLUME2->nnz_d;

						m  = malloc(NvnS[0]*Nvar * sizeof *m);  // free
						n  = malloc(NvnS[0]*Nvar * sizeof *n);  // free
						vv = malloc(NvnS[0]*Nvar * sizeof *vv); // free

						for (j = 0, jMax = NvnS[0]*Nvar; j < jMax; j++) {
							m[j]  = IndA[0]+j;
							n[j]  = IndA[0]+i;
							vv[j] = cimag(RHS_c[j])/h;
						}

						MatSetValues(*A,nnz_d,m,1,n,vv,ADD_VALUES);
						free(m); free(n); free(vv);
					} else if (assemble_type == 3) {
						if (FACET->Boundary)
							continue;

						if (side == 0) {
							if (VOLUME->indexg != FACET->VIn->indexg)
								continue;
							VOLUME2 = FACET->VOut;
							RHS_c = FACET->RHSOut_c;
						} else {
							if (VOLUME->indexg != FACET->VOut->indexg)
								continue;
							VOLUME2 = FACET->VIn;
							RHS_c = FACET->RHSIn_c;
						}

						IndA[1] = VOLUME2->IndA;
						NvnS[1] = VOLUME2->NvnS;
						nnz_d   = VOLUME2->nnz_d;

						m  = malloc(NvnS[1]*Nvar * sizeof *m);  // free
						n  = malloc(NvnS[1]*Nvar * sizeof *n);  // free
						vv = malloc(NvnS[1]*Nvar * sizeof *vv); // free

						for (j = 0, jMax = NvnS[1]*Nvar; j < jMax; j++) {
							m[j]  = IndA[1]+j;
							n[j]  = IndA[0]+i;
							vv[j] = cimag(RHS_c[j])/h;
						}

						MatSetValues(*A,nnz_d,m,1,n,vv,ADD_VALUES);
						free(m); free(n); free(vv);
					}
				}}
			}
			if (strstr(TestCase,"Poisson")) {
				VOLUME->uhat_c[i] -= h*I;
			} else {
				VOLUME->What_c[i] -= h*I;
			}
		}
	}
}

static void flag_nonzero_LHS(unsigned int *A)
{
	// Initialize DB Parameters
	unsigned int Nvar = DB.Nvar,
	             dof  = DB.dof;

	// Standard datatypes
	unsigned int i, j, iMax, jMax, side, NvnS[2], IndA[2], m, Indn;

	struct S_VOLUME *VOLUME, *VOLUME2;
	struct S_FACET  *FACET;

	for (VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next) {
		IndA[0] = VOLUME->IndA;
		NvnS[0] = VOLUME->NvnS;
		for (i = 0, iMax = NvnS[0]*Nvar; i < iMax; i++) {

			// Diagonal entries
			for (VOLUME2 = DB.VOLUME; VOLUME2; VOLUME2 = VOLUME2->next) {
				if (VOLUME->indexg != VOLUME2->indexg)
					continue;

				IndA[1] = VOLUME2->IndA;
				NvnS[1] = VOLUME2->NvnS;

				Indn = (IndA[0]+i)*dof;
				for (j = 0, jMax = NvnS[1]*Nvar; j < jMax; j++) {
					m = IndA[1]+j;
					A[Indn+m] = 1;
				}

			}

			// Off-diagonal entries
			for (FACET = DB.FACET; FACET; FACET = FACET->next) {
			for (side = 0; side < 2; side++) {
				if (FACET->Boundary)
					continue;

				if (side == 0) {
					if (VOLUME->indexg != FACET->VIn->indexg)
						continue;
					VOLUME2 = FACET->VOut;
				} else {
					if (VOLUME->indexg != FACET->VOut->indexg)
						continue;
					VOLUME2 = FACET->VIn;
				}

				IndA[1] = VOLUME2->IndA;
				NvnS[1] = VOLUME2->NvnS;

				Indn = (IndA[0]+i)*dof;
				for (j = 0, jMax = NvnS[1]*Nvar; j < jMax; j++) {
					m = IndA[1]+j;
					A[Indn+m] = 1;
				}
			}}
		}
	}
}

static void compute_A_cs_complete(Mat *A, Vec *b, Vec *x)
{
	// Initialize DB Parameters
	unsigned int Nvar = DB.Nvar,
	             dof  = DB.dof;

	// Standard datatypes
	unsigned int   i, j, iMax, jMax, NvnS[2], IndA[2], nnz_d, *A_nz;
	double         h;
	double complex *RHS_c;

	struct S_VOLUME *VOLUME, *VOLUME2;

	PetscInt    *m, *n;
	PetscScalar *vv;

	h = EPS*EPS;

	A_nz = calloc(dof*dof , sizeof *A_nz); // free

	if (*A == NULL) {
		flag_nonzero_LHS(A_nz);
		initialize_KSP(A,b,x);
	}

	for (VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next) {
		IndA[0] = VOLUME->IndA;
		NvnS[0] = VOLUME->NvnS;
		for (i = 0, iMax = NvnS[0]*Nvar; i < iMax; i++) {
			VOLUME->What_c[i] += h*I;

			explicit_VOLUME_info_c();
			explicit_FACET_info_c();
			finalize_RHS_c();

			for (VOLUME2 = DB.VOLUME; VOLUME2; VOLUME2 = VOLUME2->next) {
				IndA[1] = VOLUME2->IndA;
				if (!A_nz[(IndA[0]+i)*dof+IndA[1]])
					continue;

				NvnS[1] = VOLUME2->NvnS;
				nnz_d   = VOLUME2->nnz_d;

				m  = malloc(NvnS[1]*Nvar * sizeof *m);  // free
				n  = malloc(NvnS[1]*Nvar * sizeof *n);  // free
				vv = malloc(NvnS[1]*Nvar * sizeof *vv); // free

				RHS_c = VOLUME2->RHS_c;
				for (j = 0, jMax = NvnS[1]*Nvar; j < jMax; j++) {
					m[j]  = IndA[1]+j;
					n[j]  = IndA[0]+i;
					vv[j] = cimag(RHS_c[j])/h;
				}

				MatSetValues(*A,nnz_d,m,1,n,vv,ADD_VALUES);
				free(m); free(n); free(vv);
			}
			VOLUME->What_c[i] -= h*I;
		}
	}
	free(A_nz);

	finalize_Mat(A,1);
}

void test_integration_linearization(int nargc, char **argv)
{
	unsigned int pass;
	char         **argvNew;

	argvNew    = malloc(2          * sizeof *argvNew);  // free
	argvNew[0] = malloc(STRLEN_MAX * sizeof **argvNew); // free
	argvNew[1] = malloc(STRLEN_MAX * sizeof **argvNew); // free

	strcpy(argvNew[0],argv[0]);

	/*
	 *	Input:
	 *
	 *		Mesh file on which to run the test.
	 *
	 *	Expected Output:
	 *
	 *		The analytical linearization should match that computed using complex step.
	 *
	 */

	Mat A = NULL, A_cs = NULL, A_csc = NULL;
	Vec b = NULL, b_cs = NULL, b_csc = NULL,
	    x = NULL, x_cs = NULL, x_csc = NULL;

	// **************************************************************************************************** //
	// 2D (Mixed TRI/QUAD mesh)
	strcpy(argvNew[1],"test/Test_linearization_mixed2D");

	code_startup(nargc,argvNew,2,0);

	implicit_VOLUME_info();
	implicit_FACET_info();
//	finalize_LHS(&A,&b,&x,1);
//	finalize_LHS(&A,&b,&x,2);
//	finalize_LHS(&A,&b,&x,3);
//	finalize_Mat(&A,1);

//	compute_A_cs(&A_cs,&b_cs,&x_cs,1);
//	compute_A_cs(&A_cs,&b_cs,&x_cs,2);
//	compute_A_cs(&A_cs,&b_cs,&x_cs,3);
//	finalize_Mat(&A_cs,1);

//	MatView(A,PETSC_VIEWER_STDOUT_SELF);
//	MatView(A_cs,PETSC_VIEWER_STDOUT_SELF);
//	EXIT_MSG;

	finalize_LHS(&A,&b,&x,0);
	compute_A_cs(&A_cs,&b_cs,&x_cs,0);
	compute_A_cs_complete(&A_csc,&b_csc,&x_csc);

//	MatView(A_csc,PETSC_VIEWER_STDOUT_SELF);

	pass = 0;
	if (PetscMatAIJ_norm_diff_d(DB.dof,A,A_cs,"Inf")  < EPS &&
	    PetscMatAIJ_norm_diff_d(DB.dof,A,A_csc,"Inf") < EPS)
		pass = 1, TestDB.Npass++;

	//     0         10        20        30        40        50
	printf("Linearization (2D - Mixed):                      ");
	test_print(pass);

	finalize_ksp(&A,&b,&x,2);
	finalize_ksp(&A_cs,&b_cs,&x_cs,2);
	finalize_ksp(&A_csc,&b_csc,&x_csc,2);
	code_cleanup();

	// **************************************************************************************************** //
	// 3D (Mixed TET/PYR mesh)
	strcpy(argvNew[1],"test/Test_linearization_mixed3D_TP");

	code_startup(nargc,argvNew,2,0);

	implicit_VOLUME_info();
	implicit_FACET_info();

	finalize_LHS(&A,&b,&x,0);
	compute_A_cs(&A_cs,&b_cs,&x_cs,0);
	compute_A_cs_complete(&A_csc,&b_csc,&x_csc);

	pass = 0;
	if (PetscMatAIJ_norm_diff_d(DB.dof,A,A_cs,"Inf")  < EPS &&
	    PetscMatAIJ_norm_diff_d(DB.dof,A,A_csc,"Inf") < EPS)
		pass = 1, TestDB.Npass++;

	//     0         10        20        30        40        50
	printf("Linearization (3D - Mixed SI/PYR):               ");
	test_print(pass);

	finalize_ksp(&A,&b,&x,2);
	finalize_ksp(&A_cs,&b_cs,&x_cs,2);
	finalize_ksp(&A_csc,&b_csc,&x_csc,2);
	code_cleanup();


	// **************************************************************************************************** //
	// 3D (Mixed HEX/WEDGE mesh)
	strcpy(argvNew[1],"test/Test_linearization_mixed3D_HW");

	code_startup(nargc,argvNew,2,0);

	implicit_VOLUME_info();
	implicit_FACET_info();

	finalize_LHS(&A,&b,&x,0);
	compute_A_cs(&A_cs,&b_cs,&x_cs,0);
	compute_A_cs_complete(&A_csc,&b_csc,&x_csc);

	pass = 0;
	if (PetscMatAIJ_norm_diff_d(DB.dof,A,A_cs,"Inf")  < EPS &&
	    PetscMatAIJ_norm_diff_d(DB.dof,A,A_csc,"Inf") < EPS)
		pass = 1, TestDB.Npass++;

	//     0         10        20        30        40        50
	printf("Linearization (3D - Mixed HEX/WEDGE):            ");
	test_print(pass);

	finalize_ksp(&A,&b,&x,2);
	finalize_ksp(&A_cs,&b_cs,&x_cs,2);
	finalize_ksp(&A_csc,&b_csc,&x_csc,2);
	code_cleanup();

	free(argvNew[0]); free(argvNew[1]); free(argvNew);
}
