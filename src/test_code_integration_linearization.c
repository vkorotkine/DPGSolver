// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)
/// \file

#include "test_code_integration_linearization.h"

#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <complex.h>
#include <time.h>

#include "petscmat.h"

#include "Parameters.h"
#include "Macros.h"
#include "S_DB.h"
#include "S_VOLUME.h"
#include "S_FACE.h"
#include "Test.h"

#include "test_code_integration.h"
#include "test_support.h"

#include "solver_symmetric_functions.h" // ToBeDeleted

#include "solver.h"
#include "solver_c.h"
#include "explicit_info_c.h"

#include "explicit_VOLUME_info_c.h" // ToBeDeleted
#include "explicit_FACE_info_c.h" // ToBeDeleted
#include "finalize_RHS_c.h" // ToBeDeleted

#include "array_norm.h"

#define CHECK_VOLUME_DIAG      1
#define CHECK_VOLUME_FACE_DIAG 2
#define CHECK_VOLUME_FACE_ALL  3

/// \todo There are things to be deleted.
static void set_test_linearization_data (struct Test_Linearization *const data, char const *const TestName)
{
	data->print_mat_to_file        = false;
	data->print_timings            = false;
	data->check_full_linearization = true;
	data->check_weak_gradients     = false;

	data->PG_add        = 1;
	data->IntOrder_mult = 2;
	data->IntOrder_add  = 0;

	data->PGlobal = 3;
	data->ML      = 0;

	data->Nref        = 2;
	data->update_argv = true;

	strcpy(data->argvNew[1],TestName);
	if (strstr(TestName,"Advection")) {
		; // Do nothing
	} else if (strstr(TestName,"Poisson")) {
		; // Do nothing
	} else if (strstr(TestName,"Euler")) {
//		data->print_timings = true;
		data->update_argv = false;
	} else if (strstr(TestName,"NavierStokes")) {
//		data->print_timings = true;
		data->check_weak_gradients = true;
	} else {
		printf("%s\n",TestName); EXIT_UNSUPPORTED;
	}

	data->check_level = CHECK_VOLUME_FACE_ALL;
	if (!data->check_full_linearization) {
		data->check_level = CHECK_VOLUME_FACE_ALL;
		printf("check_level set to %d for partial contribution to the linearization.\n",data->check_level);
	}
}


static void update_VOLUME_FACEs (void);
static void perturb_solution    (void);

static void print_times         (double const *const times, bool const print_timings);
static void check_passing       (struct Test_Linearization const *const data, bool*const pass);

typedef void (*compute_A_Qhat_tdef)    (Mat A, const unsigned int check_level, const unsigned int dim);
static Mat  set_A_Qhat                 (const unsigned int check_level, const unsigned int dim,
                                        compute_A_Qhat_tdef compute_A_Qhat);
static void compute_A_Qhat             (Mat A, const unsigned int check_level, const unsigned int dim);
static void compute_A_Qhat_cs          (Mat A, const unsigned int check_level, const unsigned int dim);
static void compute_A_Qhat_cs_complete (Mat A, const unsigned int check_level, const unsigned int dim);



static void compute_A_cs               (Mat *A, const unsigned int check_level, bool const AllowOffDiag);
static void compute_A_cs_complete      (Mat *A);

static void update_times (clock_t *const tc, clock_t *const tp, double *const times, const unsigned int counter,
                          bool const print_timings);
static Mat  construct_A (void);
static void compute_A (Mat A, const unsigned int check_level);
static void assemble_A (Mat A);
static void destruct_A (Mat A);


void test_linearization
	(struct Test_Linearization *const data,
	 char const *const TestName
	)
{
	set_test_linearization_data(data,TestName);

	TestDB.PGlobal       = data->PGlobal;
	TestDB.ML            = data->ML;
	TestDB.PG_add        = data->PG_add;
	TestDB.IntOrder_add  = data->IntOrder_add;
	TestDB.IntOrder_mult = data->IntOrder_mult;

	code_startup(data->nargc,(const char*const*const) data->argvNew,data->Nref,data->update_argv);
	update_VOLUME_FACEs();
	compute_dof();

	perturb_solution();

/// \todo Move this to VOLUME struct initialization
	for (struct S_VOLUME* VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next)
		set_element(VOLUME,DB.ELEMENT);

	struct Simulation simulation = constructor_Simulation(DB.d,DB.Nvar);
	struct Context_Solver_c context_solver_c = constructor_Context_Solver_c(&simulation,DB.VOLUME);

	// Weak Gradient Linearization
	if (data->check_weak_gradients) {
		if (!strstr(TestName,"NavierStokes"))
			EXIT_UNSUPPORTED;

		data->omit_root = false;
		set_PrintName("linearization (weak gradient)",data->PrintName,&data->omit_root);
		data->omit_root = false;

		struct S_solver_info solver_info = constructor_solver_info(false,false,false,'I',DB.Method);
		compute_GradW(&solver_info,'A'); // free

		const unsigned int d = DB.d;
		bool pass = 1;
		for (size_t dim = 0; dim < d; dim++) {
			data->A     = NULL;
			data->A_cs  = NULL;
			data->A_csc = NULL;

			if (!data->check_full_linearization) {
				data->A    = set_A_Qhat(data->check_level,dim,compute_A_Qhat);
				data->A_cs = set_A_Qhat(data->check_level,dim,compute_A_Qhat_cs);
			} else {
				clock_t tp = clock(), tc;
				double times[3];
				unsigned int counter = 0;

				data->A = set_A_Qhat(data->check_level,dim,compute_A_Qhat);
				update_times(&tc,&tp,times,counter++,data->print_timings);

				data->A_cs = set_A_Qhat(data->check_level,dim,compute_A_Qhat_cs);
				update_times(&tc,&tp,times,counter++,data->print_timings);

				data->A_csc = set_A_Qhat(data->check_level,dim,compute_A_Qhat_cs_complete);
				update_times(&tc,&tp,times,counter++,data->print_timings);

				print_times(times,data->print_timings);
			}

			check_passing(data,&pass);
		}
		compute_GradW(&solver_info,'F');

		test_print2(pass,data->PrintName);
	}

	// Standard Linearization
	set_PrintName("linearization",data->PrintName,&data->omit_root);

// Move into set_A
	data->A     = construct_A(),
	data->A_cs  = construct_A(),
	data->A_csc = construct_A();

	if (!data->check_full_linearization) {
		// Note: Poisson fails symmetric for check_level = 1 and this is as expected. There is a FACE contribution from
		//       Qhat which has not yet been balanced by the FACE contribution from the 2nd equation.
		const unsigned int check_level = data->check_level;
		compute_A(data->A,check_level);
		assemble_A(data->A);

		bool AllowOffDiag = 0;
		if (check_level == CHECK_VOLUME_FACE_ALL)
			AllowOffDiag = 1;

		for (size_t i = 1; i <= check_level; i++)
			compute_A_cs(&data->A_cs,i,AllowOffDiag);
		assemble_A(data->A_cs);
	} else {
		clock_t tp = clock(), tc;
		double times[3];
		unsigned int counter = 0;

		compute_A(data->A,CHECK_VOLUME_FACE_ALL);
		assemble_A(data->A);
		update_times(&tc,&tp,times,counter++,data->print_timings);

		compute_A_cs(&data->A_cs,0,true);
		update_times(&tc,&tp,times,counter++,data->print_timings);

		compute_A_cs_complete(&data->A_csc);
		update_times(&tc,&tp,times,counter++,data->print_timings);

		print_times(times,data->print_timings);
	}

	if (data->print_mat_to_file) {
		// See commented MatView options in solver_implicit to output in more readable format.
		MatView(data->A,PETSC_VIEWER_STDOUT_SELF);
		MatView(data->A_cs,PETSC_VIEWER_STDOUT_SELF);
	}

	bool pass = 1;
	check_passing(data,&pass);
	test_print2(pass,data->PrintName);

	destructor_Simulation(&simulation);
	destructor_Context_Solver_c(&context_solver_c);

	code_cleanup();
}



// static functions

static void update_VOLUME_FACEs (void)
{
	/*
	 *	Purpose:
	 *		For each VOLUME, set the list of pointers to all adjacent FACEs.
	 */

	// Ensure that the lists are reset (in case the mesh has been updated)
	for (struct S_VOLUME *VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next) {
		for (size_t i = 0; i < NFMAX*NSUBFMAX; i++)
			VOLUME->FACE[i] = NULL;
	}

	for (struct S_FACE *FACE = DB.FACE; FACE; FACE = FACE->next) {
		update_data_FACE(FACE);

		struct S_SIDE_DATA data = FACE->data[0];
		data.VOLUME->FACE[data.Indsfh] = FACE;

		if (!FACE->Boundary) {
			data = FACE->data[1];
			data.VOLUME->FACE[data.Indsfh] = FACE;
		}
	}
}

static void perturb_solution (void)
{
	for (struct S_VOLUME *VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next) {
		const unsigned int Nvar = DB.Nvar,
						   NvnS = VOLUME->NvnS;
		for (size_t i = 0, iMax = NvnS*Nvar; i < iMax; i++)
			VOLUME->What[i] += 1e3*EPS*((double) rand() / ((double) RAND_MAX+1));
	}
}

static Mat set_A_Qhat (const unsigned int check_level, const unsigned int dim, compute_A_Qhat_tdef compute_A_Qhat)
{
	Mat A = construct_A();
	compute_A_Qhat(A,check_level,dim);
	assemble_A(A);

	return A;
}

static void check_passing(struct Test_Linearization const *const data, bool*const pass)
{
	double diff_cs = 0.0;
	if (data->check_full_linearization)
		diff_cs = PetscMatAIJ_norm_diff_d(DB.dof,data->A_cs,data->A_csc,"Inf",false);

	if (diff_cs > 2e1*EPS)
		*pass = 0;

	double const diff_LHS = PetscMatAIJ_norm_diff_d(DB.dof,data->A_cs,data->A,"Inf",true);
	if (diff_LHS > 2e1*EPS)
		*pass = 0;

	if (DB.Symmetric) {
		PetscBool Symmetric = 0;
		MatIsSymmetric(data->A,1e1*EPS,&Symmetric);
		if (!Symmetric) {
			*pass = 0;
			printf("Failed symmetric.\n");
		}
	}

	if (!(*pass)) {
		if (data->check_full_linearization)
			printf("diff_cs:  %e\n",diff_cs);
		printf("diff_LHS: %e\n",diff_LHS);
	}

	destruct_A(data->A);
	destruct_A(data->A_cs);
	destruct_A(data->A_csc);
}

static void update_times (clock_t *const tc, clock_t *const tp, double *const times, const unsigned int counter,
                          bool const print_timings)
{
	if (!print_timings)
		return;

	*tc = clock();
	times[counter] = (*tc-*tp)/(1.0*CLOCKS_PER_SEC);
	*tp = *tc;
}

static void print_times (double const *const times, bool const print_timings)
{
	if (!print_timings)
		return;

	printf("Timing ratios: ");
	for (size_t i = 0; i < 2; i++)
		printf("% .2f ",times[i+1]/times[0]);
	printf("\n");
}

static Mat construct_A (void)
{
	struct S_solver_info solver_info = constructor_solver_info(false,false,false,'I',DB.Method);
	solver_info.create_RHS = false;
	initialize_petsc_structs(&solver_info);

	return solver_info.A;
}

static void compute_A (Mat A, const unsigned int check_level)
{
	struct S_solver_info solver_info = constructor_solver_info(false,false,false,'I',DB.Method);
	solver_info.create_RHS = false;
	solver_info.A = A;

	if (check_level == CHECK_VOLUME_DIAG)
		solver_info.compute_F = false;

	compute_RLHS(&solver_info);
}

static void assemble_A (Mat A)
{
	struct S_solver_info solver_info = constructor_solver_info(false,false,false,'I',DB.Method);
	solver_info.create_RHS = false;
	solver_info.A = A;
	assemble_petsc_structs(&solver_info);
}

static void destruct_A (Mat A)
{
	struct S_solver_info solver_info = constructor_solver_info(false,false,false,'I',DB.Method);
	solver_info.create_RHS = false;
	solver_info.A = A;
	destroy_petsc_structs(&solver_info);
}

static void flag_nonzero_LHS(unsigned int *A)
{
	const unsigned int Nvar = DB.Nvar,
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
					if (VOLUME->indexg != FACE->VL->indexg)
						continue;
					VOLUME2 = FACE->VR;
				} else {
					if (VOLUME->indexg != FACE->VR->indexg)
						continue;
					VOLUME2 = FACE->VL;
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

static void compute_A_cs(Mat *const A, const unsigned int check_level, bool const AllowOffDiag)
{
	bool check_off_diagonal = false;
	if (AllowOffDiag) {
		if (DB.Viscous)
			check_off_diagonal = true;
	}

	if (check_level == 0) {
		compute_A_cs(A,1,AllowOffDiag);
		compute_A_cs(A,2,AllowOffDiag);
		compute_A_cs(A,3,AllowOffDiag);

		assemble_A(*A);
		return;
	}

	const unsigned int Nvar = DB.Nvar;
	unsigned int       NvnS[2], IndA[2];
	for (struct S_VOLUME *VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next) {
		NvnS[0] = VOLUME->NvnS;

		// Note: Initialize with zeros for linear cases.
		if (strstr(DB.TestCase,"Advection") || strstr(DB.TestCase,"Poisson")) {
			if (VOLUME->What_c)
				free(VOLUME->What_c);

			VOLUME->What_c = calloc(NvnS[0]*Nvar , sizeof *(VOLUME->What_c));
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
			VOLUME->What_c[i] += CX_STEP*I;

			struct S_solver_info solver_info = constructor_solver_info(false,false,false,'E',DB.Method);
			solver_info.VOLUME_perturbed = VOLUME;

			compute_GradW_c(&solver_info,'A'); // free
			if (check_level == 1) {
				explicit_VOLUME_info_c(VOLUME,0);
				correct_collocated_for_symmetry_c(VOLUME,0,true,false);
			} else if (check_level == 2 || check_level == 3) {
				explicit_FACE_info_c(VOLUME,0);
				correct_collocated_for_symmetry_c(VOLUME,false,false,true);
			} else {
				EXIT_UNSUPPORTED;
			}
			compute_GradW_c(&solver_info,'F');

			if (check_level == 1) {
				const unsigned int dof = DB.dof;

				unsigned int *A_nz  = NULL;
				if (check_off_diagonal) {
					A_nz = calloc(dof*dof , sizeof *A_nz); // free
					flag_nonzero_LHS(A_nz);
				}

				for (struct S_VOLUME *VOLUME2 = DB.VOLUME; VOLUME2; VOLUME2 = VOLUME2->next) {
					if (!check_off_diagonal) {
						if (VOLUME->indexg != VOLUME2->indexg)
							continue;
					}

					IndA[1] = VOLUME2->IndA;
					if (check_off_diagonal && !A_nz[(IndA[0]+i)*dof+IndA[1]])
						continue;

					NvnS[1] = VOLUME2->NvnS;
					const unsigned int nnz_d = VOLUME2->nnz_d;

					PetscInt    *const m  = malloc(NvnS[1]*Nvar * sizeof *m);  // free
					PetscInt    *const n  = malloc(NvnS[1]*Nvar * sizeof *n);  // free
					PetscScalar *const vv = malloc(NvnS[1]*Nvar * sizeof *vv); // free

					double complex const *const RHS_c = VOLUME2->RHS_c;
					for (size_t j = 0, jMax = NvnS[1]*Nvar; j < jMax; j++) {
						m[j]  = IndA[1]+j;
						n[j]  = IndA[0]+i;
						vv[j] = cimag(RHS_c[j])/CX_STEP;
					}

					MatSetValues(*A,nnz_d,m,1,n,vv,ADD_VALUES);
					free(m); free(n); free(vv);
				}

				if (check_off_diagonal)
					free(A_nz);
			} else if (check_level == 2 || check_level == 3) {
				for (struct S_FACE *FACE = DB.FACE; FACE; FACE = FACE->next) {
				for (size_t side = 0; side < 2; side++) {
					double complex const *RHS_c;

					struct S_VOLUME *VOLUME2;
					if (check_level == 2) {
						if (side == 0) {
							VOLUME2 = FACE->VL;
							RHS_c = FACE->RHSL_c;
						} else {
							if (FACE->Boundary)
								continue;
							VOLUME2 = FACE->VR;
							RHS_c = FACE->RHSR_c;
						}
						if (VOLUME->indexg != VOLUME2->indexg)
							continue;

						IndA[1] = VOLUME2->IndA;
						NvnS[1] = VOLUME2->NvnS;
						const unsigned int nnz_d = VOLUME2->nnz_d;

						PetscInt    *const m  = malloc(NvnS[0]*Nvar * sizeof *m);  // free
						PetscInt    *const n  = malloc(NvnS[0]*Nvar * sizeof *n);  // free
						PetscScalar *const vv = malloc(NvnS[0]*Nvar * sizeof *vv); // free

						for (size_t j = 0, jMax = NvnS[0]*Nvar; j < jMax; j++) {
							m[j]  = IndA[0]+j;
							n[j]  = IndA[0]+i;
							vv[j] = cimag(RHS_c[j])/CX_STEP;
						}

						MatSetValues(*A,nnz_d,m,1,n,vv,ADD_VALUES);
						free(m); free(n); free(vv);
					} else if (check_level == 3) {
						if (FACE->Boundary)
							continue;

						if (side == 0 && VOLUME == FACE->VL) {
							VOLUME2 = FACE->VR;
							RHS_c = FACE->RHSR_c;
						} else if (side == 1 && VOLUME == FACE->VR) {
							VOLUME2 = FACE->VL;
							RHS_c = FACE->RHSL_c;
						} else {
							continue;
						}

						IndA[1] = VOLUME2->IndA;
						NvnS[1] = VOLUME2->NvnS;
						const unsigned int nnz_d = VOLUME2->nnz_d;

						PetscInt    *const m  = malloc(NvnS[1]*Nvar * sizeof *m);  // free
						PetscInt    *const n  = malloc(NvnS[1]*Nvar * sizeof *n);  // free
						PetscScalar *const vv = malloc(NvnS[1]*Nvar * sizeof *vv); // free

						for (size_t j = 0, jMax = NvnS[1]*Nvar; j < jMax; j++) {
							m[j]  = IndA[1]+j;
							n[j]  = IndA[0]+i;
							vv[j] = cimag(RHS_c[j])/CX_STEP;
						}

						MatSetValues(*A,nnz_d,m,1,n,vv,ADD_VALUES);
						free(m); free(n); free(vv);
					}
				}}
			} else {
				EXIT_UNSUPPORTED;
			}
			VOLUME->What_c[i] -= CX_STEP*I;
		}
	}
}

static void compute_A_cs_complete(Mat *A)
{
	const unsigned int dof = DB.dof;

	unsigned int *A_nz = calloc(dof*dof , sizeof *A_nz); // free
	flag_nonzero_LHS(A_nz);

	for (struct S_VOLUME *VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next) {
		const unsigned int Nvar = DB.Nvar;
		unsigned int       NvnS[2], IndA[2];

		IndA[0] = VOLUME->IndA;
		NvnS[0] = VOLUME->NvnS;
		for (size_t i = 0, iMax = NvnS[0]*Nvar; i < iMax; i++) {
			VOLUME->What_c[i] += CX_STEP*I;

			struct S_solver_info solver_info = constructor_solver_info(false,false,false,'E',DB.Method);
			solver_info.VOLUME_perturbed = VOLUME;

			compute_GradW_c(&solver_info,'A'); // free
			explicit_VOLUME_info_c(VOLUME,0);
			explicit_FACE_info_c(VOLUME,0);
			compute_GradW_c(&solver_info,'F');
			finalize_RHS_c(VOLUME,0);

			// correct_FACE = false as FACE terms were added to VOLUME RHS in finalize_RHS_c
			correct_collocated_for_symmetry_c(VOLUME,0,true,false);

			for (struct S_VOLUME *VOLUME2 = DB.VOLUME; VOLUME2; VOLUME2 = VOLUME2->next) {
				IndA[1] = VOLUME2->IndA;
				if (!A_nz[(IndA[0]+i)*dof+IndA[1]])
					continue;

				NvnS[1] = VOLUME2->NvnS;
				const unsigned int nnz_d = VOLUME2->nnz_d;

				PetscInt    *const m  = malloc(NvnS[1]*Nvar * sizeof *m);  // free
				PetscInt    *const n  = malloc(NvnS[1]*Nvar * sizeof *n);  // free
				PetscScalar *const vv = malloc(NvnS[1]*Nvar * sizeof *vv); // free

				double complex const *const RHS_c = VOLUME2->RHS_c;
				for (size_t j = 0, jMax = NvnS[1]*Nvar; j < jMax; j++) {
					m[j]  = IndA[1]+j;
					n[j]  = IndA[0]+i;
					vv[j] = cimag(RHS_c[j])/CX_STEP;
				}

				MatSetValues(*A,nnz_d,m,1,n,vv,ADD_VALUES);
				free(m); free(n); free(vv);
			}

			VOLUME->What_c[i] -= CX_STEP*I;
		}
	}
	free(A_nz);

	assemble_A(*A);
}

static void fill_PetscMat_GradW (const double*const LHS, const struct S_VOLUME*const V0, const struct S_VOLUME*const V1,
                                 const bool boundary, Mat A)

{
	/*
	 *	Purpose:
	 *		Add individual Qhat_What contributions to the Petsc Mat.
	 *
	 *	Comments:
	 *		This function is nearly identical to the standard fill_PetscMat function but allows for different indexing
	 *		of the LHS array due to the redundant of Qhat*_What terms when not on a boundary. These terms being
	 *		redundant resulted in memory not being allocated for them in the Qhat*_What arrays.
	 */

	struct S_LHS_info LHS_info = constructor_LHS_info((double*) LHS,V0,V1,ADD_VALUES,false);

	const size_t*const IndA = LHS_info.IndA,
	            *const Nn   = LHS_info.Nn;

	PetscInt idxm[Nn[0]],
	         idxn[Nn[1]];
	const InsertMode addv = LHS_info.addv;

	const unsigned int Neq  = DB.Neq,
	                   Nvar = DB.Nvar;

	const double*const zeros = calloc(Nn[0]*Nn[1] , sizeof *zeros); // free
	for (size_t eq = 0; eq < Neq; eq++) {
		const size_t Indm = IndA[0] + eq*Nn[0];
		for (size_t i = 0; i < Nn[0]; i++)
			idxm[i] = Indm+i;

		for (size_t var = 0; var < Nvar; var++) {
			const size_t Indn = IndA[1] + var*Nn[1];
			for (size_t i = 0; i < Nn[1]; i++)
				idxn[i] = Indn+i;

			const PetscScalar* vv = NULL;
			if (boundary)
				vv = &LHS[(eq*Nvar+var)*Nn[0]*Nn[1]];
			else
				vv = ( eq == var ? &LHS[0] : zeros );

			MatSetValues(A,Nn[0],idxm,Nn[1],idxn,vv,addv);
		}
	}
	free((double*) zeros);
}

static void fill_PetscMat_GradW_c (const double complex*const RHS_c, const struct S_VOLUME*const V0,
                                   const struct S_VOLUME*const V1, const size_t IndA_col, Mat A)
{
	/*
	 *	Purpose:
	 *		Add individual Qhat_c contributions to the Petsc Mat.
	 *
	 *	Comments:
	 *		This function builds A column-wise using the complex Qhat terms.
	 */

	const unsigned int Nvar = DB.Nvar;

	const unsigned int IndA[] = { V0->IndA, V1->IndA },
	                   Nn[]   = { V0->NvnS*Nvar, V1->NvnS*Nvar };

	PetscInt idxm[Nn[0]];
	PetscInt idxn[] = { IndA_col };

	PetscScalar*const vv = malloc(Nn[0] * sizeof *vv); // free
	for (size_t i = 0; i < Nn[0]; i++) {
		idxm[i] = IndA[0]+i;
		vv[i]   = cimag(RHS_c[i])/CX_STEP;
	}

	MatSetValues(A,Nn[0],idxm,1,idxn,vv,ADD_VALUES);
	free(vv);
}

static void compute_A_Qhat(Mat A, const unsigned int check_level, const unsigned int dim)
{
	/*
	 *	Purpose:
	 *		Assembles the various contributions to the weak gradient in A.
	 *
	 *	Comments:
	 *		check_level is a flag for which parts of the global matrix are checked.
	 */

	switch (check_level) {
	case 3: // off-diagonal contributions
		for (struct S_FACE *FACE = DB.FACE; FACE; FACE = FACE->next) {
			if (FACE->Boundary)
				continue;

			fill_PetscMat_GradW(FACE->QhatL_WhatR[dim],FACE->VL,FACE->VR,false,A);
			fill_PetscMat_GradW(FACE->QhatR_WhatL[dim],FACE->VR,FACE->VL,false,A);
		}
		// fallthrough
	case 2: // diagonal FACE contributions
		for (struct S_FACE *FACE = DB.FACE; FACE; FACE = FACE->next) {
			fill_PetscMat_GradW(FACE->QhatL_WhatL[dim],FACE->VL,FACE->VL,FACE->Boundary,A);

			if (FACE->Boundary)
				continue;

			fill_PetscMat_GradW(FACE->QhatR_WhatR[dim],FACE->VR,FACE->VR,false,A);
		}
		// fallthrough
	case 1: // diagonal VOLUME contributions
		for (struct S_VOLUME *VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next) {
			fill_PetscMat_GradW(VOLUME->QhatV_What[dim],VOLUME,VOLUME,false,A);
		}
		break;
	default:
		EXIT_UNSUPPORTED;
		break;
	}
}

static void initialize_What_c (void)
{
	const unsigned int Nvar = DB.Nvar;
	for (struct S_VOLUME *VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next) {
		const unsigned int NvnS = VOLUME->NvnS;

// ToBeModified (Should previously have been freed elsewhere?)
		if (VOLUME->What_c)
			free(VOLUME->What_c);
		VOLUME->What_c = malloc(NvnS*Nvar * sizeof *(VOLUME->What_c));

		for (size_t i = 0, iMax = NvnS*Nvar; i < iMax; i++)
			VOLUME->What_c[i] = VOLUME->What[i];
	}
}

static void compute_A_Qhat_cs(Mat A, const unsigned int check_level, const unsigned int dim)
{
	initialize_What_c();

	const unsigned int Nvar = DB.Nvar;
	for (struct S_VOLUME *VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next) {
		const unsigned int IndA = VOLUME->IndA,
		                   Nn   = VOLUME->NvnS*Nvar;
		for (size_t i = 0; i < Nn; i++) {
			VOLUME->What_c[i] += CX_STEP*I;

			struct S_solver_info solver_info = constructor_solver_info(false,false,false,'E',DB.Method);
			solver_info.VOLUME_perturbed = VOLUME;

			compute_GradW_c(&solver_info,'A'); // free

			switch (check_level) {
			case 3:
				for (struct S_FACE *FACE = DB.FACE; FACE; FACE = FACE->next) {
					if (FACE->Boundary)
						continue;

					if (VOLUME == FACE->VL)
						fill_PetscMat_GradW_c(FACE->QhatR_c[dim],FACE->VR,FACE->VL,IndA+i,A);
					if (VOLUME == FACE->VR)
						fill_PetscMat_GradW_c(FACE->QhatL_c[dim],FACE->VL,FACE->VR,IndA+i,A);
				}
				// fallthrough
			case 2:
				for (struct S_FACE *FACE = DB.FACE; FACE; FACE = FACE->next) {
					if (VOLUME == FACE->VL)
						fill_PetscMat_GradW_c(FACE->QhatL_c[dim],FACE->VL,FACE->VL,IndA+i,A);

					if (FACE->Boundary)
						continue;

					if (VOLUME == FACE->VR)
						fill_PetscMat_GradW_c(FACE->QhatR_c[dim],FACE->VR,FACE->VR,IndA+i,A);
				}
				// fallthrough
			case 1:
				for (struct S_VOLUME *VOLUME2 = DB.VOLUME; VOLUME2; VOLUME2 = VOLUME2->next) {
					if (VOLUME != VOLUME2)
						continue;

					fill_PetscMat_GradW_c(VOLUME2->QhatV_c[dim],VOLUME2,VOLUME,IndA+i,A);
				}
				break;
			default:
				EXIT_UNSUPPORTED;
				break;
			}
			compute_GradW_c(&solver_info,'F');

			VOLUME->What_c[i] -= CX_STEP*I;
		}
	}
}

static void compute_A_Qhat_cs_complete (Mat A, const unsigned int check_level, const unsigned int dim)
{
	if (check_level != CHECK_VOLUME_FACE_ALL)
		EXIT_UNSUPPORTED;

	const unsigned int Nvar = DB.Nvar;
	for (struct S_VOLUME *VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next) {
		const unsigned int IndA = VOLUME->IndA,
		                   Nn   = VOLUME->NvnS*Nvar;
		for (size_t i = 0; i < Nn; i++) {
			VOLUME->What_c[i] += CX_STEP*I;

			struct S_solver_info solver_info = constructor_solver_info(false,false,false,'E',DB.Method);
			solver_info.VOLUME_perturbed = VOLUME;

			compute_GradW_c(&solver_info,'A'); // free
			compute_GradW_c(&solver_info,'F');

			struct S_LOCAL_MESH_ELEMENTS local_ELEMENTs = compute_local_ELEMENT_list(VOLUME,'V');

			for (struct S_VOLUME *VOLUME2 = DB.VOLUME; VOLUME2; VOLUME2 = VOLUME2->next) {
				if (!is_VOLUME_in_local_list(VOLUME2,&local_ELEMENTs))
					continue;

				fill_PetscMat_GradW_c(VOLUME2->Qhat_c[dim],VOLUME2,VOLUME,IndA+i,A);
			}

			VOLUME->What_c[i] -= CX_STEP*I;
		}
	}
}
