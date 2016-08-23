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
#include "explicit_VOLUME_info_c.h"
#include "implicit_VOLUME_info.h"
#include "explicit_FACET_info_c.h"
#include "implicit_FACET_info.h"
#include "finalize_LHS.h"

/*
 *	Purpose:
 *		Test correctness of implementation of code linearization.
 *
 *	Comments:
 *
 *	Notation:
 *
 *	References:
 */

static void compute_A_cs(Mat *A, Vec *b, const unsigned int assemble_type)
{
	// Initialize DB Paramters
	unsigned int Nvar = DB.Nvar;

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
		initialize_KSP(A,b);

	for (VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next) {
		NvnS[0] = VOLUME->NvnS;

		if (VOLUME->What_c)
			free(VOLUME->What_c);

		VOLUME->What_c = malloc(NvnS[0]*Nvar * sizeof *(VOLUME->What_c));

		for (i = 0, iMax = NvnS[0]*Nvar; i < iMax; i++)
			VOLUME->What_c[i] = VOLUME->What[i];
	}

	for (VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next) {
		IndA[0] = VOLUME->IndA;
		NvnS[0] = VOLUME->NvnS;
		for (i = 0, iMax = NvnS[0]*Nvar; i < iMax; i++) {
			VOLUME->What_c[i] += h*I;

			switch (assemble_type) {
			default: // 0
				explicit_VOLUME_info_c();
				explicit_FACET_info_c();
				printf("Error: Not yet implemented.\n"), EXIT_MSG;
				break;
			case 1:
				explicit_VOLUME_info_c();
				break;
			case 2:
			case 3:
				explicit_FACET_info_c();
				printf("Error: Not yet implemented.\n"), EXIT_MSG;
				break;
			}

			if (assemble_type == 0 || assemble_type == 1) {
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

					MatSetValues(*A,nnz_d,m,1,n,vv,INSERT_VALUES);

					free(m); free(n); free(vv);
				}
			}

			if (assemble_type == 0 || assemble_type == 2 || assemble_type == 3) {
				for (FACET = DB.FACET; FACET; FACET = FACET->next) {
				for (side = 0; side < 2; side++) {
					if (side == 0)
						VOLUME2 = FACET->VIn;
					else
						VOLUME2 = FACET->VOut;

					if (VOLUME->indexg == VOLUME2->indexg) {
						IndA[1] = VOLUME2->IndA;
						NvnS[1] = VOLUME2->NvnS;
						nnz_d   = VOLUME2->nnz_d;

						m  = malloc(NvnS[1]*Nvar * sizeof *m);  // free
						n  = malloc(NvnS[1]*Nvar * sizeof *n);  // free
						vv = malloc(NvnS[1]*Nvar * sizeof *vv); // free

						if (assemble_type == 0 || assemble_type == 2) {
							if (side == 0)
								RHS_c = FACET->RHSIn_c;
							else
								RHS_c = FACET->RHSOut_c;

							for (j = 0, jMax = NvnS[1]*Nvar; j < jMax; j++) {
								m[j]  = IndA[1]+j;
								n[j]  = IndA[0]+i;
								vv[j] = cimag(RHS_c[j])/h;
							}

							MatSetValues(*A,nnz_d,m,1,n,vv,INSERT_VALUES);
						}

						if (assemble_type == 0 || assemble_type == 3) {
							if (side == 0)
								RHS_c = FACET->RHSOut_c;
							else
								RHS_c = FACET->RHSIn_c;

							for (j = 0, jMax = NvnS[1]*Nvar; j < jMax; j++) {
								m[j]  = IndA[1]+j;
								n[j]  = IndA[0]+i;
								vv[j] = cimag(RHS_c[j])/h;
							}

							MatSetValues(*A,nnz_d,m,1,n,vv,INSERT_VALUES);
						}

						free(m); free(n); free(vv);
					}
				}}
			}
			VOLUME->What_c[i] -= h*I;
		}
	}
	MatAssemblyBegin(*A,MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(*A,MAT_FINAL_ASSEMBLY);
//	MatView(*A,PETSC_VIEWER_STDOUT_SELF);
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

	Mat A = NULL, A_cs = NULL;
	Vec b = NULL, b_cs = NULL;

	// **************************************************************************************************** //
	// LINEs


	// **************************************************************************************************** //
	// TRIs
	strcpy(argvNew[1],"test/Test_linearization_TRI");

	code_startup(nargc,argvNew,2);

	implicit_VOLUME_info();

implicit_FACET_info();
finalize_LHS(&A,&b,1);
compute_A_cs(&A_cs,&b_cs,2);

	finalize_LHS(&A,&b,1);

	compute_A_cs(&A_cs,&b_cs,1);

	pass = 0;
	if (PetscMatAIJ_norm_diff_d(DB.dof,A,A_cs,"Inf") < EPS)
		pass = 1, TestDB.Npass++;
	//     0         10        20        30        40        50
	printf("Linearization - VOLUME diag (TRI  ):             ");
	test_print(pass);


	// Don't forget to MatDestroy b when implemented (ToBeDeleted)
	MatDestroy(&A);    A    = NULL;
	MatDestroy(&A_cs); A_cs = NULL;

	code_cleanup(0);



	free(argvNew[0]); free(argvNew[1]); free(argvNew);
}
