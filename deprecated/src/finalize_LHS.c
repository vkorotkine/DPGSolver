// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#include "finalize_LHS.h"

#include <stdlib.h>
#include <stdio.h>

#include "petscvec.h"
#include "petscmat.h"

#include "Macros.h"
#include "S_DB.h"
#include "S_VOLUME.h"
#include "S_FACE.h"

#include "finalize_RHS.h"
#include "solver_symmetric_functions.h"
#include "array_print.h"

/*
 *	Purpose:
 *		Finalize LHS term by summing VOLUME and FACE contributions.
 *
 *	Comments:
 *
 *	Notation:
 *
 *	References:
 */

void initialize_KSP(Mat *A, Vec *b, Vec *x)
{
	/*
	 *	Comments:
	 *		Potentially need to add in matrix for preconditioner.
	 *		Likely generalize this to use MPI matrix once MPI is pursued.
	 *		Potentially need to add CHKERRQ for MatCreate. (ToBeDeleted)
	 */
EXIT_UNSUPPORTED;
printf("iKSP\n");

	unsigned int const dof = DB.dof;

	MPI_Comm comm = MPI_COMM_WORLD;
	PetscInt *nnz = malloc(dof * sizeof *nnz); // free

	size_t Indi = 0;
	for (struct S_VOLUME *VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next) {
		unsigned int const nnz_d = VOLUME->nnz_d;

		for (size_t i = 0; i < nnz_d; i++)
			nnz[Indi+i] = nnz_d + VOLUME->nnz_o;

		Indi += nnz_d;
	}

	MatCreateSeqAIJ(comm,dof,dof,0,nnz,A); // keep (requires external MatDestroy)
	VecCreateSeq(comm,dof,b);              // keep (requires external VecDestroy)
	VecCreateSeq(comm,dof,x);              // keep (requires external VecDestroy)

	free(nnz);
}

void finalize_Vec(Vec *a, const unsigned int finalize_type)
{
	switch (finalize_type) {
		case 1: // Petsc assemble
			VecAssemblyBegin(*a);
			VecAssemblyEnd(*a);
			break;
		case 2: // Petsc destroy
			VecDestroy(a);
			a = NULL;
			break;
		default:
			printf("Error: Unsupported finalize_type.\n"), EXIT_MSG;
			break;
	}
}

void finalize_Mat(Mat *A, const unsigned int finalize_type)
{
	switch (finalize_type) {
		case 1: // Petsc assemble
			MatAssemblyBegin(*A,MAT_FINAL_ASSEMBLY);
			MatAssemblyEnd(*A,MAT_FINAL_ASSEMBLY);
			break;
		case 2: // Petsc destroy
			MatDestroy(A);
			A = NULL;
			break;
		default:
			EXIT_UNSUPPORTED;
			break;
	}
}

void finalize_ksp(Mat *A, Vec *b, Vec *x, const unsigned int finalize_type)
{
	finalize_Mat(A,finalize_type);
	finalize_Vec(b,finalize_type);
	finalize_Vec(x,finalize_type);
}

void assemble_RHS(Vec *b, Vec *x)
{
	// Initialize DB Parameters
	unsigned int Nvar = DB.Nvar;

	// Standard datatypes
	unsigned int i, iMax, Indb, NvnS;
	double       *RHS;

	struct S_VOLUME *VOLUME;

	PetscInt    *ix;
	PetscScalar *y;

	for (VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next) {
		Indb = VOLUME->IndA;
		NvnS = VOLUME->NvnS;
		RHS  = VOLUME->RHS;
array_print_d(NvnS,Nvar,RHS,'C');

		iMax = NvnS*Nvar;
		ix = malloc(iMax * sizeof *ix); // free
		y  = malloc(iMax * sizeof *y);  // free

		for (i = 0; i < iMax; i++) {
			ix[i] = Indb+i;
			y[i]  = -(*RHS++);
		}

		VecSetValues(*b,iMax,ix,y,INSERT_VALUES);
		VecSetValues(*x,iMax,ix,y,INSERT_VALUES);

		free(ix);
		free(y);
	}
}

static void compute2_dof(void);
double finalize_LHS(Mat *A, Vec *b, Vec *x, const unsigned int assemble_type)
{
	/*
	 *	Comments:
	 *		assemble_type is a flag for which parts of the global matrix are assembled.
	 */

	compute2_dof();

	if (*A == NULL)
		initialize_KSP(A,b,x);

	double maxRHS = 1.0;
	switch (assemble_type) {
	default: // 0
		maxRHS = finalize_RHS();
		correct_collocated_for_symmetry();

		assemble_RHS(b,x);

		finalize_ksp(A,b,x,1);
		break;
	case 1: // diagonal VOLUME contributions
	case 2: // diagonal FACE contributions
	case 3: // off-diagonal contributions
		// Do nothing
		break;
	}
	return maxRHS;
}

static void compute2_dof(void)
{
	/*
	 *	Comments:
	 *		Likely put in a reordering algorithm based on conclusions of Persson (2008). (ToBeDeleted)
	 *
	 *	References:
	 *		Persson(2008)-Newton-GMRES_Preconditioning_for_Discontinuous_Galerkin_Discretizations_of_the_Navier-Stokes_
	 *		              Equations
	 */

	unsigned int Nvar = DB.Nvar;

	unsigned int dof = 0;
	for (struct S_VOLUME *VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next) {
		VOLUME->IndA = dof;

		unsigned int nnz_d = Nvar*(VOLUME->NvnS);
		dof += nnz_d;

		VOLUME->nnz_d = nnz_d;
		VOLUME->nnz_o = 0;
	}

	for (struct S_FACE *FACE = DB.FACE; FACE; FACE = FACE->next) {
		struct S_VOLUME *VL = FACE->VL,
		                *VR = FACE->VR;

		if (VL->indexg != VR->indexg) {
			VL->nnz_o += VR->nnz_d;
			VR->nnz_o += VL->nnz_d;
		}
	}

	DB.dof = dof;
}
