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
#include "array_print.h"

/*
 *	Purpose:
 *		Provide functions for linearization integration testing.
 *
 *	Comments:
 *		The linearization is verified by comparing with the output when using the complex step method.
 *
 *		For second order equations, the verification of the linearization of the weak gradients is also performed.
 *
 *		For second order equations, the VOLUME->RHS_c contributes off-diagonal terms to the global system matrix due to
 *		the use of the fully corrected gradient. A flag is provided to avoid the computation of these terms when
 *		checking only the diagonal VOLUME/FACE contributions to the global system using compute_A_cs with assembly_type
 *		equal to 1 or 2. The current implementation for the Poisson equation avoids this problem by separating the
 *		contributions of the local and correction contributions to the VOLUME gradient and this issue is thus not
 *		present in that case.
 *
 *		A similar formulation for the Navier-Stokes equations may be possible due to the linearity of the viscous flux
 *		with respect to the solution gradients. INVESTIGATE. (ToBeModified)
 *
 *	Notation:
 *
 *	References:
 *		Martins(2003)-The_Complex-Step_Derivative_Approximation
 */

static void compute_A_cs               (Mat *A, Vec *b, Vec *x, unsigned int const assemble_type,
                                        bool const AllowOffDiag);
static void compute_A_cs_complete      (Mat *A, Vec *b, Vec *x);
static void finalize_LHS_Qhat          (Mat *const A, Vec *const b, Vec *const x, unsigned int const assemble_type,
                                        unsigned int const dim);
static void compute_A_Qhat_cs          (Mat *const A, Vec *const b, Vec *const x, unsigned int const assemble_type,
                                        unsigned int const dim);
static void compute_A_Qhat_cs_complete (Mat *const A, Vec *const b, Vec *const x, unsigned int const dim);

static void set_test_linearization_data(struct S_linearization *const data, char const *const TestName)
{
	// default values
	TestDB.CheckOffDiagonal = 0;

	data->PrintEnabled           = 0;
	data->CheckFullLinearization = 1;
	data->CheckWeakGradients     = 0;
	data->CheckSymmetric         = 0;

	data->PG_add        = 1;
	data->IntOrder_mult = 2;
	data->IntOrder_add  = 0;

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
		data->CheckWeakGradients = 1;
		if (strstr(TestName,"TRI")) {
			strcpy(data->argvNew[1],"test/NavierStokes/Test_NavierStokes_TaylorCouette_ToBeCurvedTRI");
		} else {
			EXIT_UNSUPPORTED;
		}
	} else {
		printf("%s\n",TestName); EXIT_UNSUPPORTED;
	}
}


static void check_passing(struct S_linearization const *const data, unsigned int *pass)
{
	double diff_cs = 0.0;
	if (data->CheckFullLinearization)
		diff_cs = PetscMatAIJ_norm_diff_d(DB.dof,data->A_cs,data->A_csc,"Inf");

	if (diff_cs > 1e1*EPS)
		*pass = 0;

	double const diff_LHS = PetscMatAIJ_norm_diff_d(DB.dof,data->A_cs,data->A,"Inf");
	if (diff_LHS > 1e1*EPS)
		*pass = 0;

	if (data->CheckSymmetric) {
		PetscBool Symmetric = 0;
		MatIsSymmetric(data->A,1e2*EPS,&Symmetric);
//		MatIsSymmetric(data->A_cs,EPS,&Symmetric);
		if (!Symmetric) {
			*pass = 0;
			printf("Failed symmetric.\n");
		}
	}

	if (!(*pass)) {
		if (data->CheckFullLinearization)
			printf("diff_cs:  %e\n",diff_cs);
		printf("diff_LHS: %e\n",diff_LHS);
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

	bool const CheckFullLinearization = data->CheckFullLinearization,
	           CheckWeakGradients     = data->CheckWeakGradients,
	           PrintEnabled           = data->PrintEnabled;

	unsigned int const Nref        = data->Nref,
	                   update_argv = data->update_argv;

	TestDB.PGlobal       = data->PGlobal;
	TestDB.ML            = data->ML;
	TestDB.PG_add        = data->PG_add;
	TestDB.IntOrder_add  = data->IntOrder_add;
	TestDB.IntOrder_mult = data->IntOrder_mult;

	code_startup(nargc,argvNew,Nref,update_argv);

	// Perturb the solution for nonlinear equations
	if (strstr(TestName,"Euler") || strstr(TestName,"NavierStokes")) {
		for (struct S_VOLUME *VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next) {
			unsigned int const Nvar = DB.Nvar,
			                   NvnS = VOLUME->NvnS;
			for (size_t i = 0, iMax = NvnS*Nvar; i < iMax; i++)
				VOLUME->What[i] += 1e3*EPS*((double) rand() / ((double) RAND_MAX+1));
		}
	} else if (strstr(TestName,"Poisson")) {
		; // Do nothing
	} else {
		EXIT_UNSUPPORTED;
	}

	for (size_t nTest = 0; nTest < 2; nTest++) {
		if (nTest == 0) { // Check weak gradients
			if (!CheckWeakGradients)
				continue;

			unsigned int const d = DB.d;

			Mat A[DMAX] = { NULL }, A_cs[DMAX] = { NULL }, A_csc[DMAX] = { NULL };
			Vec b[DMAX] = { NULL }, b_cs[DMAX] = { NULL }, b_csc[DMAX] = { NULL },
			    x[DMAX] = { NULL }, x_cs[DMAX] = { NULL }, x_csc[DMAX] = { NULL };

			set_PrintName("linearization (weak gradient)",data->PrintName,&data->TestTRI);
			if (strstr(TestName,"NavierStokes")) {
				implicit_GradW();
			} else {
				EXIT_UNSUPPORTED;
			}

			if (!CheckFullLinearization) {
				for (size_t dim = 0; dim < d; dim++) {
					unsigned int CheckLevel = 3;
					for (size_t i = 1; i <= CheckLevel; i++)
						finalize_LHS_Qhat(&A[dim],&b[dim],&x[dim],i,dim);
					finalize_Mat(&A[dim],1);

					for (size_t i = 1; i <= CheckLevel; i++)
						compute_A_Qhat_cs(&A_cs[dim],&b_cs[dim],&x_cs[dim],i,dim);
					finalize_Mat(&A_cs[dim],1);
				}
			} else {
				for (size_t dim = 0; dim < d; dim++) {
					finalize_LHS_Qhat(&A[dim],&b[dim],&x[dim],0,dim);
					compute_A_Qhat_cs(&A_cs[dim],&b_cs[dim],&x_cs[dim],0,dim);
					compute_A_Qhat_cs_complete(&A_csc[dim],&b_csc[dim],&x_csc[dim],dim);
				}
			}

			unsigned int pass = 1;
			for (size_t dim = 0; dim < d; dim++) {
				data->A = A[dim]; data->A_cs = A_cs[dim]; data->A_csc = A_csc[dim];
				data->b = b[dim]; data->b_cs = b_cs[dim]; data->b_csc = b_csc[dim];
				data->x = x[dim]; data->x_cs = x_cs[dim]; data->x_csc = x_csc[dim];

				check_passing(data,&pass);
			}
			test_print2(pass,data->PrintName);

			for (size_t dim = 0; dim < d; dim++) {
				finalize_ksp(&A[dim],&b[dim],&x[dim],2);
				finalize_ksp(&A_cs[dim],&b_cs[dim],&x_cs[dim],2);
				finalize_ksp(&A_csc[dim],&b_csc[dim],&x_csc[dim],2);
			}
		} else { // Standard (Full linearization)
			Mat A = NULL, A_cs = NULL, A_csc = NULL;
			Vec b = NULL, b_cs = NULL, b_csc = NULL, x = NULL, x_cs = NULL, x_csc = NULL;

			set_PrintName("linearization",data->PrintName,&data->TestTRI);
			if (strstr(TestName,"Poisson")) {
				implicit_info_Poisson();
			} else if (strstr(TestName,"Euler") || strstr(TestName,"NavierStokes")) {
				// Comment implicit_FACE_info here and used Q_vI = ChiS_vI*VOLUME->QhatV to only check VOLUME-VOLUME
				// contribution to the VOLUME term.

				implicit_GradW(); // Only executed if DB.Viscous = 1
				implicit_VOLUME_info();
				implicit_FACE_info();
			} else {
				EXIT_UNSUPPORTED;
			}

			if (!CheckFullLinearization) {
				unsigned int const CheckLevel = 3;
				for (size_t i = 1; i <= CheckLevel; i++)
					finalize_LHS(&A,&b,&x,i);
				finalize_Mat(&A,1);

				bool AllowOffDiag = 0;
				if (CheckLevel == 3)
					AllowOffDiag = 1;

				for (size_t i = 1; i <= CheckLevel; i++)
					compute_A_cs(&A_cs,&b_cs,&x_cs,i,AllowOffDiag);
				finalize_Mat(&A_cs,1);
			} else {
				finalize_LHS(&A,&b,&x,0);
				compute_A_cs(&A_cs,&b_cs,&x_cs,0,1);
				compute_A_cs_complete(&A_csc,&b_csc,&x_csc);
			}

			if (PrintEnabled) {
				MatView(A,PETSC_VIEWER_STDOUT_SELF);
				MatView(A_cs,PETSC_VIEWER_STDOUT_SELF);
			}

			data->A = A; data->A_cs = A_cs; data->A_csc = A_csc;
			data->b = b; data->b_cs = b_cs; data->b_csc = b_csc;
			data->x = x; data->x_cs = x_cs; data->x_csc = x_csc;

			unsigned int pass = 1;
			check_passing(data,&pass);
			test_print2(pass,data->PrintName);

			finalize_ksp(&A,&b,&x,2);
			finalize_ksp(&A_cs,&b_cs,&x_cs,2);
			finalize_ksp(&A_csc,&b_csc,&x_csc,2);
		}
	}
	code_cleanup();
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

static void compute_A_cs(Mat *const A, Vec *const b, Vec *const x, unsigned int const assemble_type,
                         bool const AllowOffDiag)
{
	if (*A == NULL)
		initialize_KSP(A,b,x);

	if (AllowOffDiag) {
		if (strstr(DB.TestCase,"NavierStokes")) {
			TestDB.CheckOffDiagonal = 1;
		} else if (strstr(DB.TestCase,"Poisson") || strstr(DB.TestCase,"Euler")) {
			// Note: If modified to be more similar to the Navier-Stokes solver, the Poisson solver would also likely
			//       need to used the treatment for the Navier-Stokes below (i.e. there would be off-diagonal terms in
			//       the global linear system matrix when only running considering the VOLUME terms (due to the use of
			//       the fully corrected gradient). As it is currently implemented, these contributions appear to be
			//       included in FACE->LHS(LL/RR), meaning that this is not necessary.
			; // Do nothing
		} else {
			EXIT_UNSUPPORTED;
		}
	}

	if (assemble_type == 0) {
		compute_A_cs(A,b,x,1,AllowOffDiag);
		compute_A_cs(A,b,x,2,AllowOffDiag);
		compute_A_cs(A,b,x,3,AllowOffDiag);

		finalize_Mat(A,1);
		return;
	}

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
				explicit_VOLUME_info_c();
				explicit_FACE_info_c();
			} else {
				EXIT_UNSUPPORTED;
			}

			if (assemble_type == 1) {
				unsigned int const dof = DB.dof;

				bool CheckOffDiagonal = TestDB.CheckOffDiagonal;
				unsigned int *A_nz  = NULL;
				if (CheckOffDiagonal) {
					A_nz = calloc(dof*dof , sizeof *A_nz); // free
					flag_nonzero_LHS(A_nz);
				}

				for (struct S_VOLUME *VOLUME2 = DB.VOLUME; VOLUME2; VOLUME2 = VOLUME2->next) {
					if (!CheckOffDiagonal) {
						if (VOLUME->indexg != VOLUME2->indexg)
							continue;
					}

					IndA[1] = VOLUME2->IndA;
					if (CheckOffDiagonal && !A_nz[(IndA[0]+i)*dof+IndA[1]])
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

				if (CheckOffDiagonal)
					free(A_nz);
			} else if (assemble_type == 2 || assemble_type == 3) {
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

						if (side == 0 && VOLUME == FACE->VIn) {
							VOLUME2 = FACE->VOut;
							RHS_c = FACE->RHSOut_c;
						} else if (side == 1 && VOLUME == FACE->VOut) {
							VOLUME2 = FACE->VIn;
							RHS_c = FACE->RHSIn_c;
						} else {
							continue;
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

static void finalize_LHS_Qhat(Mat *const A, Vec *const b, Vec *const x, unsigned int const assemble_type,
                              unsigned int const dim)
{
	/*
	 *	Purpose:
	 *		Assembles the various contributions to the weak gradient in A.
	 *
	 *	Comments:
	 *		assemble_type is a flag for which parts of the global matrix are assembled.
	 */

	unsigned int Nvar = DB.Nvar,
	             Neq  = DB.Neq;

	compute_dof();

	if (*A == NULL)
		initialize_KSP(A,b,x);

	switch (assemble_type) {
	default: // 0
		finalize_LHS_Qhat(A,b,x,1,dim);
		finalize_LHS_Qhat(A,b,x,2,dim);
		finalize_LHS_Qhat(A,b,x,3,dim);

		finalize_Mat(A,1);
		break;
	case 1: // diagonal VOLUME contributions
		for (struct S_VOLUME *VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next) {
			unsigned int const IndA = VOLUME->IndA,
			                   NvnS = VOLUME->NvnS;

			PetscInt *const m = malloc(NvnS * sizeof *m), // free
			         *const n = malloc(NvnS * sizeof *n); // free

			double *zeros = calloc(NvnS*NvnS , sizeof *zeros); // free

			for (size_t eq = 0; eq < Neq; eq++) {
				size_t const Indm = IndA + eq*NvnS;
				for (size_t i = 0; i < NvnS; i++)
					m[i] = Indm+i;

				for (size_t var = 0; var < Nvar; var++) {
					size_t const Indn = IndA + var*NvnS;
					for (size_t i = 0; i < NvnS; i++)
						n[i] = Indn+i;

					PetscScalar const *vv;
					if (eq == var)
						vv = VOLUME->QhatV_What[dim];
					else
						vv = zeros;

					MatSetValues(*A,NvnS,m,NvnS,n,vv,ADD_VALUES);
				}
			}
			free(m);
			free(n);
			free(zeros);
		}
		break;
	case 2: // diagonal FACE contributions
		for (struct S_FACE *FACE = DB.FACE; FACE; FACE = FACE->next) {
		for (size_t side = 0; side < 2; side++) {
			double const *QhatF_What = NULL;

			struct S_VOLUME const *VOLUME;
			if (side == 0) {
				VOLUME = FACE->VIn;
				QhatF_What = FACE->Qhat_WhatLL[dim];
			} else {
				if (FACE->Boundary)
					continue;
				VOLUME = FACE->VOut;
				QhatF_What = FACE->Qhat_WhatRR[dim];
			}

			unsigned int IndA = VOLUME->IndA,
			             NvnS = VOLUME->NvnS;

			PetscInt *const m = malloc(NvnS * sizeof *m), // free
			         *const n = malloc(NvnS * sizeof *n); // free

			double *zeros = calloc(NvnS*NvnS , sizeof *zeros); // free

			for (size_t eq = 0; eq < Neq; eq++) {
				size_t const Indm = IndA + eq*NvnS;
				for (size_t i = 0; i < NvnS; i++)
					m[i] = Indm+i;

				for (size_t var = 0; var < Nvar; var++) {
					size_t const Indn = IndA + var*NvnS;
					for (size_t i = 0; i < NvnS; i++)
						n[i] = Indn+i;

					PetscScalar const *vv;
					if (FACE->Boundary) {
						vv = &QhatF_What[(eq*Nvar+var)*NvnS*NvnS];
					} else {
						if (eq == var)
							vv = &QhatF_What[0];
						else
							vv = zeros;
					}

					MatSetValues(*A,NvnS,m,NvnS,n,vv,ADD_VALUES);
				}
			}
			free(m);
			free(n);
			free(zeros);
		}}
		break;
	case 3: // off-diagonal contributions
		for (struct S_FACE *FACE = DB.FACE; FACE; FACE = FACE->next) {
		for (size_t side = 0; side < 2; side++) {
			if (FACE->Boundary)
				continue;

			double const *QhatF_What = NULL;

			struct S_VOLUME const *VOLUME, *VOLUME2;
			if (side == 0) {
				VOLUME  = FACE->VOut;
				VOLUME2 = FACE->VIn;
				QhatF_What = FACE->Qhat_WhatLR[dim];
			} else {
				VOLUME  = FACE->VIn;
				VOLUME2 = FACE->VOut;
				QhatF_What = FACE->Qhat_WhatRL[dim];
			}

			unsigned int const IndA  = VOLUME->IndA,
			                   IndA2 = VOLUME2->IndA,
			                   NvnS  = VOLUME->NvnS,
			                   NvnS2 = VOLUME2->NvnS;

			PetscInt *const m = malloc(NvnS  * sizeof *m), // free
			         *const n = malloc(NvnS2 * sizeof *n); // free

			double *zeros = calloc(NvnS*NvnS2 , sizeof *zeros); // free

			for (size_t eq = 0; eq < Neq; eq++) {
				size_t const Indm = IndA + eq*NvnS;
				for (size_t i = 0; i < NvnS; i++)
					m[i] = Indm+i;

				for (size_t var = 0; var < Nvar; var++) {
					size_t const Indn = IndA2 + var*NvnS2;
					for (size_t i = 0; i < NvnS2; i++)
						n[i] = Indn+i;

					PetscScalar const *vv;
					if (eq == var)
						vv = &QhatF_What[0];
					else
						vv = zeros;

					MatSetValues(*A,NvnS,m,NvnS2,n,vv,ADD_VALUES);
				}
			}
			free(m); free(n);
		}}
		break;
	}
}

static void compute_A_Qhat_cs(Mat *const A, Vec *const b, Vec *const x, unsigned int const assemble_type,
                              unsigned int const dim)
{
	if (*A == NULL)
		initialize_KSP(A,b,x);

	if (assemble_type == 0) {
		compute_A_Qhat_cs(A,b,x,1,dim);
		compute_A_Qhat_cs(A,b,x,2,dim);
		compute_A_Qhat_cs(A,b,x,3,dim);

		finalize_Mat(A,1);
		return;
	}

	unsigned int const Nvar = DB.Nvar;
	unsigned int       NvnS[2], IndA[2];
	for (struct S_VOLUME *VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next) {
		NvnS[0] = VOLUME->NvnS;

		// Note: Initialize with zeros for linear cases.
		if (strstr(DB.TestCase,"NavierStokes")) {
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
			if (strstr(DB.TestCase,"NavierStokes")) {
				VOLUME->What_c[i] += h*I;

				explicit_GradW_c();
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

					double complex const *const QhatV_c = VOLUME2->QhatV_c[dim];
					for (size_t j = 0, jMax = NvnS[1]*Nvar; j < jMax; j++) {
						m[j]  = IndA[1]+j;
						n[j]  = IndA[0]+i;
						vv[j] = cimag(QhatV_c[j])/h;
					}

					MatSetValues(*A,nnz_d,m,1,n,vv,ADD_VALUES);
					free(m); free(n); free(vv);
				}
			}

			if (assemble_type == 2 || assemble_type == 3) {
				for (struct S_FACE *FACE = DB.FACE; FACE; FACE = FACE->next) {
				for (size_t side = 0; side < 2; side++) {
					double complex const *QhatF_c;

					struct S_VOLUME *VOLUME2;
					if (assemble_type == 2) {
						if (side == 0) {
							VOLUME2 = FACE->VIn;
							QhatF_c = FACE->QhatL_c[dim];
						} else {
							if (FACE->Boundary)
								continue;
							VOLUME2 = FACE->VOut;
							QhatF_c = FACE->QhatR_c[dim];
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
							vv[j] = cimag(QhatF_c[j])/h;
						}

						MatSetValues(*A,nnz_d,m,1,n,vv,ADD_VALUES);
						free(m); free(n); free(vv);
					} else if (assemble_type == 3) {
						if (FACE->Boundary)
							continue;

						if (side == 0 && VOLUME == FACE->VIn) {
							VOLUME2 = FACE->VOut;
							QhatF_c = FACE->QhatR_c[dim];
						} else if (side == 1 && VOLUME == FACE->VOut) {
							VOLUME2 = FACE->VIn;
							QhatF_c = FACE->QhatL_c[dim];
						} else {
							continue;
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
							vv[j] = cimag(QhatF_c[j])/h;
						}

						MatSetValues(*A,nnz_d,m,1,n,vv,ADD_VALUES);
						free(m); free(n); free(vv);
					}
				}}
			}
			if (strstr(DB.TestCase,"NavierStokes")) {
				VOLUME->What_c[i] -= h*I;
			} else {
				EXIT_UNSUPPORTED;
			}
		}
	}
}

static void compute_A_Qhat_cs_complete(Mat *const A, Vec *const b, Vec *const x, unsigned int const dim)
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
			if (strstr(DB.TestCase,"NavierStokes")) {
				VOLUME->What_c[i] += h*I;

				explicit_GradW_c();
			} else {
				EXIT_UNSUPPORTED;
			}

			for (struct S_VOLUME *VOLUME2 = DB.VOLUME; VOLUME2; VOLUME2 = VOLUME2->next) {
				IndA[1] = VOLUME2->IndA;
				if (!A_nz[(IndA[0]+i)*dof+IndA[1]])
					continue;

				NvnS[1] = VOLUME2->NvnS;
				unsigned int const nnz_d = VOLUME2->nnz_d;

				PetscInt    *const m  = malloc(NvnS[1]*Nvar * sizeof *m);  // free
				PetscInt    *const n  = malloc(NvnS[1]*Nvar * sizeof *n);  // free
				PetscScalar *const vv = malloc(NvnS[1]*Nvar * sizeof *vv); // free

				double complex const *const Qhat_c = VOLUME2->Qhat_c[dim];
				for (size_t j = 0, jMax = NvnS[1]*Nvar; j < jMax; j++) {
					m[j]  = IndA[1]+j;
					n[j]  = IndA[0]+i;
					vv[j] = cimag(Qhat_c[j])/h;
				}

				MatSetValues(*A,nnz_d,m,1,n,vv,ADD_VALUES);
				free(m); free(n); free(vv);
			}

			if (strstr(DB.TestCase,"NavierStokes")) {
				VOLUME->What_c[i] -= h*I;
			} else {
				EXIT_UNSUPPORTED;
			}
		}
	}
	free(A_nz);

	finalize_Mat(A,1);
}
