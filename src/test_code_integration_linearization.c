// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#include "test_code_integration_linearization.h"

#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <complex.h>

#include "petscmat.h"

#include "Parameters.h"
#include "Macros.h"
#include "S_DB.h"
#include "S_VOLUME.h"
#include "S_FACE.h"
#include "Test.h"

#include "test_code_integration.h"
#include "test_support.h"
#include "solver_Poisson.h"
#include "finalize_LHS.h"
#include "solver_Poisson_c.h"
#include "explicit_VOLUME_info_c.h"
#include "implicit_VOLUME_info.h"
#include "explicit_FACE_info_c.h"
#include "implicit_FACE_info.h"
#include "explicit_GradW_c.h"
#include "implicit_GradW.h"
#include "finalize_RHS_c.h"
#include "array_norm.h"

/*
 *	Purpose:
 *		Provide functions for linearization integration testing.
 *
 *	Comments:
 *		The linearization is verified by comparing with the output when using the complex step method.
 *
 *	Notation:
 *
 *	References:
 *		Martins(2003)-The_Complex-Step_Derivative_Approximation
 */

static void compute_A_cs          (Mat *A, Vec *b, Vec *x, unsigned int const assemble_type);
static void compute_A_cs_complete (Mat *A, Vec *b, Vec *x);

static void set_test_linearization_data(struct S_linearization *const data, char const *const TestName)
{
	// default values
	data->PrintEnabled           = 0;
	data->CheckFullLinearization = 1;
	data->CheckSymmetric         = 0;

	data->PG_add        = 1;
	data->IntOrder_mult = 2;

	data->PGlobal = 3;
	data->ML      = 0;

	data->Nref        = 2;
	data->update_argv = 1;

	if (strstr(TestName,"Poisson")) {
		data->CheckSymmetric = 1;
		if (strstr(TestName,"MIXED2D")) {
			strcpy(data->argvNew[1],"test/Poisson/Test_Poisson_n-Ball_HollowSection_CurvedMIXED2D");
		} else {
			// 3D TET test previously had Nref = 0, update_argv = 1
			EXIT_UNSUPPORTED;
		}
	} else if (strstr(TestName,"Euler")) {
		data->update_argv = 0;
// Modify this to be consisten with other PDE types (i.e. test/Euler/Test_Euler... (ToBeDeleted)
		if (strstr(TestName,"MIXED2D")) {
			strcpy(data->argvNew[1],"test/linearization/Test_linearization_ToBeCurvedMIXED2D");
		} else if (strstr(TestName,"MIXED_TET_PYR")) {
			strcpy(data->argvNew[1],"test/linearization/Test_linearization_ToBeCurvedMIXED3D_TP");
		} else if (strstr(TestName,"MIXED_HEX_WEDGE")) {
			strcpy(data->argvNew[1],"test/linearization/Test_linearization_ToBeCurvedMIXED3D_HW");
		} else {
			EXIT_UNSUPPORTED;
		}
	} else if (strstr(TestName,"NavierStokes")) {
		if (strstr(TestName,"TRI")) {
			strcpy(data->argvNew[1],"test/NavierStokes/Test_NavierStokes_TaylorCouette_ToBeCurvedTRI");
		} else {
			EXIT_UNSUPPORTED;
		}
	} else {
		printf("%s\n",TestName); EXIT_UNSUPPORTED;
	}
}

void test_linearization(struct S_linearization *const data, char const *const TestName)
{
	/*
	 *	Expected Output:
	 *		Correspondence of LHS matrices computed using complex step and exact linearization.
	 *
	 *	Comments
	 *
	 *		By default, the complete linearization is checked here. However, linearizations of the individual
	 *		contributions listed below may be checked separately:
	 *			0) complete check (default)
	 *			1) diagonal VOLUME contributions
	 *			2) diagonal FACE contributions
	 *			3) off-diagonal FACE contributions
	 *		Further, the complete linearization can be checked either through the assembly of the individual
	 *		contributions or more simply by using the assembled RHS directly.
	 */

	set_test_linearization_data(data,TestName);

	int  const               nargc   = data->nargc;
	char const *const *const argvNew = (char const *const *const) data->argvNew;

	bool const CheckSymmetric         = data->CheckSymmetric,
	           CheckFullLinearization = data->CheckFullLinearization,
	           PrintEnabled           = data->PrintEnabled;

	unsigned int const Nref        = data->Nref,
	                   update_argv = data->update_argv;

	TestDB.PGlobal       = data->PGlobal;
	TestDB.ML            = data->ML;
	TestDB.PG_add        = data->PG_add;
	TestDB.IntOrder_mult = data->IntOrder_mult;

	Mat A, A_cs, A_csc;
	Vec b, b_cs, b_csc, x, x_cs, x_csc;

	A = data->A; A_cs = data->A_cs; A_csc = data->A_csc;
	b = data->b; b_cs = data->b_cs; b_csc = data->b_csc;
	x = data->x; x_cs = data->x_cs; x_csc = data->x_csc;

	code_startup(nargc,argvNew,Nref,update_argv);

	if (strstr(TestName,"Poisson")) {
		implicit_info_Poisson();
	} else if (strstr(TestName,"Euler") || strstr(TestName,"NavierStokes")) {
		implicit_GradW(); // Only run if DB.Viscous = 1
		implicit_VOLUME_info();
		implicit_FACE_info();
	} else {
		EXIT_UNSUPPORTED;
	}

	if (!CheckFullLinearization) {
		finalize_LHS(&A,&b,&x,1);
//		finalize_LHS(&A,&b,&x,2);
//		finalize_LHS(&A,&b,&x,3);
		finalize_Mat(&A,1);

		compute_A_cs(&A_cs,&b_cs,&x_cs,1);
//		compute_A_cs(&A_cs,&b_cs,&x_cs,2);
//		compute_A_cs(&A_cs,&b_cs,&x_cs,3);
		finalize_Mat(&A_cs,1);
	} else {
		finalize_LHS(&A,&b,&x,0);
		compute_A_cs(&A_cs,&b_cs,&x_cs,0);
		compute_A_cs_complete(&A_csc,&b_csc,&x_csc);
	}

	if (PrintEnabled) {
		MatView(A,PETSC_VIEWER_STDOUT_SELF);
		MatView(A_cs,PETSC_VIEWER_STDOUT_SELF);
	}

	unsigned int pass = 0;
	if (PetscMatAIJ_norm_diff_d(DB.dof,A_cs,A,"Inf")     < 1e2*EPS &&
	    PetscMatAIJ_norm_diff_d(DB.dof,A_cs,A_csc,"Inf") < 1e2*EPS) {
		if (!CheckSymmetric)
			pass = 1;
		else {
			PetscBool Symmetric = 0;
			MatIsSymmetric(A,1e2*EPS,&Symmetric);
//			MatIsSymmetric(A_cs,EPS,&Symmetric);
			if (Symmetric)
				pass = 1;
			else
				printf("Failed symmetric.\n");
		}
	} else {
		printf("%e %e\n",PetscMatAIJ_norm_diff_d(DB.dof,A_cs,A,"Inf"),PetscMatAIJ_norm_diff_d(DB.dof,A_cs,A_csc,"Inf"));
	}

	set_PrintName("linearization",data->PrintName,&data->TestTRI);
	test_print2(pass,data->PrintName);

	finalize_ksp(&A,&b,&x,2);
	finalize_ksp(&A_cs,&b_cs,&x_cs,2);
	finalize_ksp(&A_csc,&b_csc,&x_csc,2);
	code_cleanup();
}

static void compute_A_cs(Mat *const A, Vec *const b, Vec *const x, unsigned int const assemble_type)
{
	if (!assemble_type) {
		compute_A_cs(A,b,x,1);
		compute_A_cs(A,b,x,2);
		compute_A_cs(A,b,x,3);

		finalize_Mat(A,1);
		return;
	}

	if (*A == NULL)
		initialize_KSP(A,b,x);

	unsigned int const Nvar = DB.Nvar;
	unsigned int       NvnS[2], IndA[2];
	for (struct S_VOLUME *VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next) {
		NvnS[0] = VOLUME->NvnS;

		// Note: Initialize with zeros for linear cases.
		if (strstr(DB.TestCase,"Poisson")) {
			if (VOLUME->uhat_c)
				free(VOLUME->uhat_c);

			VOLUME->uhat_c = calloc(NvnS[0]*Nvar , sizeof *(VOLUME->uhat_c));
		} else if (strstr(DB.TestCase,"Euler") || strstr(DB.TestCase,"NavierStokes")) {
			if (VOLUME->What_c)
				free(VOLUME->What_c);
			VOLUME->What_c = malloc(NvnS[0]*Nvar * sizeof *(VOLUME->What_c));

			for (size_t i = 0, iMax = NvnS[0]*Nvar; i < iMax; i++)
				VOLUME->What_c[i] = VOLUME->What[i];
		} else {
			EXIT_UNSUPPORTED;
		}
	}

	for (struct S_VOLUME *VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next) {
		IndA[0] = VOLUME->IndA;
		NvnS[0] = VOLUME->NvnS;
		for (size_t i = 0, iMax = NvnS[0]*Nvar; i < iMax; i++) {
			double const h = EPS*EPS;
			if (strstr(DB.TestCase,"Poisson")) {
				VOLUME->uhat_c[i] += h*I;

				compute_qhat_VOLUME_c();
				compute_qhat_FACE_c();
				switch (assemble_type) {
				default: // 0
					compute_uhat_VOLUME_c();
					compute_uhat_FACE_c();
					break;
				case 1:
					compute_uhat_VOLUME_c();
					break;
				case 2:
				case 3:
					compute_uhat_FACE_c();
					break;
				}
			} else if (strstr(DB.TestCase,"Euler") || strstr(DB.TestCase,"NavierStokes")) {
				VOLUME->What_c[i] += h*I;

				explicit_GradW_c();
				switch (assemble_type) {
				default: // 0
					explicit_VOLUME_info_c();
					explicit_FACE_info_c();
					break;
				case 1:
					explicit_VOLUME_info_c();
					break;
				case 2:
				case 3:
					explicit_FACE_info_c();
					break;
				}
			} else {
				EXIT_UNSUPPORTED;
			}

			if (assemble_type == 1) {
				for (struct S_VOLUME *VOLUME2 = DB.VOLUME; VOLUME2; VOLUME2 = VOLUME2->next) {
					if (VOLUME->indexg != VOLUME2->indexg)
						continue;

					IndA[1] = VOLUME2->IndA;
					NvnS[1] = VOLUME2->NvnS;
					unsigned int const nnz_d = VOLUME2->nnz_d;

					PetscInt    *const m  = malloc(NvnS[1]*Nvar * sizeof *m);  // free
					PetscInt    *const n  = malloc(NvnS[1]*Nvar * sizeof *n);  // free
					PetscScalar *const vv = malloc(NvnS[1]*Nvar * sizeof *vv); // free

					double complex const *const RHS_c = VOLUME2->RHS_c;
					for (size_t j = 0, jMax = NvnS[1]*Nvar; j < jMax; j++) {
						m[j]  = IndA[1]+j;
						n[j]  = IndA[0]+i;
						vv[j] = cimag(RHS_c[j])/h;
					}

					MatSetValues(*A,nnz_d,m,1,n,vv,ADD_VALUES);
					free(m); free(n); free(vv);
				}
			}

			if (assemble_type == 2 || assemble_type == 3) {
				for (struct S_FACE *FACE = DB.FACE; FACE; FACE = FACE->next) {
				for (size_t side = 0; side < 2; side++) {
					double complex const *RHS_c;

					struct S_VOLUME *VOLUME2;
					if (assemble_type == 2) {
						if (side == 0) {
							VOLUME2 = FACE->VIn;
							RHS_c = FACE->RHSIn_c;
						} else {
							if (FACE->Boundary)
								continue;
							VOLUME2 = FACE->VOut;
							RHS_c = FACE->RHSOut_c;
						}
						if (VOLUME->indexg != VOLUME2->indexg)
							continue;

						IndA[1] = VOLUME2->IndA;
						NvnS[1] = VOLUME2->NvnS;
						unsigned int const nnz_d = VOLUME2->nnz_d;

						PetscInt    *const m  = malloc(NvnS[0]*Nvar * sizeof *m);  // free
						PetscInt    *const n  = malloc(NvnS[0]*Nvar * sizeof *n);  // free
						PetscScalar *const vv = malloc(NvnS[0]*Nvar * sizeof *vv); // free

						for (size_t j = 0, jMax = NvnS[0]*Nvar; j < jMax; j++) {
							m[j]  = IndA[0]+j;
							n[j]  = IndA[0]+i;
							vv[j] = cimag(RHS_c[j])/h;
						}

						MatSetValues(*A,nnz_d,m,1,n,vv,ADD_VALUES);
						free(m); free(n); free(vv);
					} else if (assemble_type == 3) {
						if (FACE->Boundary)
							continue;

						if (side == 0) {
							if (VOLUME->indexg != FACE->VIn->indexg)
								continue;
							VOLUME2 = FACE->VOut;
							RHS_c = FACE->RHSOut_c;
						} else {
							if (VOLUME->indexg != FACE->VOut->indexg)
								continue;
							VOLUME2 = FACE->VIn;
							RHS_c = FACE->RHSIn_c;
						}

						IndA[1] = VOLUME2->IndA;
						NvnS[1] = VOLUME2->NvnS;
						unsigned int const nnz_d = VOLUME2->nnz_d;

						PetscInt    *const m  = malloc(NvnS[1]*Nvar * sizeof *m);  // free
						PetscInt    *const n  = malloc(NvnS[1]*Nvar * sizeof *n);  // free
						PetscScalar *const vv = malloc(NvnS[1]*Nvar * sizeof *vv); // free

						for (size_t j = 0, jMax = NvnS[1]*Nvar; j < jMax; j++) {
							m[j]  = IndA[1]+j;
							n[j]  = IndA[0]+i;
							vv[j] = cimag(RHS_c[j])/h;
						}

						MatSetValues(*A,nnz_d,m,1,n,vv,ADD_VALUES);
						free(m); free(n); free(vv);
					}
				}}
			}
			if (strstr(DB.TestCase,"Poisson")) {
				VOLUME->uhat_c[i] -= h*I;
			} else if (strstr(DB.TestCase,"Euler") || strstr(DB.TestCase,"NavierStokes")) {
				VOLUME->What_c[i] -= h*I;
			} else {
				EXIT_UNSUPPORTED;
			}
		}
	}
}

static void flag_nonzero_LHS(unsigned int *A)
{
	unsigned int const Nvar = DB.Nvar,
	                   dof  = DB.dof;
	unsigned int       NvnS[2], IndA[2];
	for (struct S_VOLUME *VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next) {
		IndA[0] = VOLUME->IndA;
		NvnS[0] = VOLUME->NvnS;
		for (size_t i = 0, iMax = NvnS[0]*Nvar; i < iMax; i++) {

			// Diagonal entries
			for (struct S_VOLUME *VOLUME2 = DB.VOLUME; VOLUME2; VOLUME2 = VOLUME2->next) {
				if (VOLUME->indexg != VOLUME2->indexg)
					continue;

				IndA[1] = VOLUME2->IndA;
				NvnS[1] = VOLUME2->NvnS;

				size_t const Indn = (IndA[0]+i)*dof;
				for (size_t j = 0, jMax = NvnS[1]*Nvar; j < jMax; j++) {
					size_t const m = IndA[1]+j;
					A[Indn+m] = 1;
				}

			}

			// Off-diagonal entries
			for (struct S_FACE *FACE = DB.FACE; FACE; FACE = FACE->next) {
			for (size_t side = 0; side < 2; side++) {
				struct S_VOLUME *VOLUME2;

				if (FACE->Boundary)
					continue;

				if (side == 0) {
					if (VOLUME->indexg != FACE->VIn->indexg)
						continue;
					VOLUME2 = FACE->VOut;
				} else {
					if (VOLUME->indexg != FACE->VOut->indexg)
						continue;
					VOLUME2 = FACE->VIn;
				}

				IndA[1] = VOLUME2->IndA;
				NvnS[1] = VOLUME2->NvnS;

				size_t const Indn = (IndA[0]+i)*dof;
				for (size_t j = 0, jMax = NvnS[1]*Nvar; j < jMax; j++) {
					size_t const m = IndA[1]+j;
					A[Indn+m] = 1;
				}
			}}
		}
	}
}

static void compute_A_cs_complete(Mat *A, Vec *b, Vec *x)
{
	unsigned int const dof = DB.dof;

	unsigned int *A_nz = calloc(dof*dof , sizeof *A_nz); // free
	if (*A == NULL) {
		flag_nonzero_LHS(A_nz);
		initialize_KSP(A,b,x);
	}

	unsigned int const Nvar = DB.Nvar;
	unsigned int       NvnS[2], IndA[2];
	for (struct S_VOLUME *VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next) {
		IndA[0] = VOLUME->IndA;
		NvnS[0] = VOLUME->NvnS;
		for (size_t i = 0, iMax = NvnS[0]*Nvar; i < iMax; i++) {
			double const h = EPS*EPS;
			if (strstr(DB.TestCase,"Poisson")) {
				VOLUME->uhat_c[i] += h*I;

				compute_qhat_VOLUME_c();
				compute_qhat_FACE_c();
				compute_uhat_VOLUME_c();
				compute_uhat_FACE_c();
			} else if (strstr(DB.TestCase,"Euler") || strstr(DB.TestCase,"NavierStokes")) {
				VOLUME->What_c[i] += h*I;

				explicit_GradW_c();
				explicit_VOLUME_info_c();
				explicit_FACE_info_c();
			} else {
				EXIT_UNSUPPORTED;
			}
			finalize_RHS_c();

			for (struct S_VOLUME *VOLUME2 = DB.VOLUME; VOLUME2; VOLUME2 = VOLUME2->next) {
				IndA[1] = VOLUME2->IndA;
				if (!A_nz[(IndA[0]+i)*dof+IndA[1]])
					continue;

				NvnS[1] = VOLUME2->NvnS;
				unsigned int const nnz_d = VOLUME2->nnz_d;

				PetscInt    *const m  = malloc(NvnS[1]*Nvar * sizeof *m);  // free
				PetscInt    *const n  = malloc(NvnS[1]*Nvar * sizeof *n);  // free
				PetscScalar *const vv = malloc(NvnS[1]*Nvar * sizeof *vv); // free

				double complex const *const RHS_c = VOLUME2->RHS_c;
				for (size_t j = 0, jMax = NvnS[1]*Nvar; j < jMax; j++) {
					m[j]  = IndA[1]+j;
					n[j]  = IndA[0]+i;
					vv[j] = cimag(RHS_c[j])/h;
				}

				MatSetValues(*A,nnz_d,m,1,n,vv,ADD_VALUES);
				free(m); free(n); free(vv);
			}

			if (strstr(DB.TestCase,"Poisson")) {
				VOLUME->uhat_c[i] -= h*I;
			} else if (strstr(DB.TestCase,"Euler") || strstr(DB.TestCase,"NavierStokes")) {
				VOLUME->What_c[i] -= h*I;
			} else {
				EXIT_UNSUPPORTED;
			}
		}
	}
	free(A_nz);

	finalize_Mat(A,1);
}
