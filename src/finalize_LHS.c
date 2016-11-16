// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

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

	// Initialize DB Parameters
	unsigned int dof = DB.dof;

	// Standard datatypes
	unsigned int i, nnz_d, Indi;

	struct S_VOLUME *VOLUME;

	MPI_Comm comm = MPI_COMM_WORLD;
	PetscInt *nnz;

	nnz = malloc(dof * sizeof *nnz); // free

	Indi = 0;
	for (VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next) {
		nnz_d = VOLUME->nnz_d;

		for (i = 0; i < nnz_d; i++)
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
			printf("Error: Unsupported finalize_type.\n"), EXIT_MSG;
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

double finalize_LHS(Mat *A, Vec *b, Vec *x, const unsigned int assemble_type)
{
	/*
	 *	Comments:
	 *		assemble_type is a flag for which parts of the global matrix are assembled.
	 */

	// Initialize DB Parameters
	unsigned int Nvar = DB.Nvar,
	             Neq  = DB.Neq;

	// Standard datatypes
	unsigned int i, NvnS, NvnS2, eq, var, side,
	             IndA, IndA2, Indm, Indn;
	double       maxRHS, *LHS;

	struct S_VOLUME *VOLUME, *VOLUME2;
	struct S_FACE  *FACE;

	PetscInt    *m, *n;
	PetscScalar *vv;

	compute_dof();

	if (*A == NULL)
		initialize_KSP(A,b,x);

	maxRHS = 1.0;
	switch (assemble_type) {
	default: // 0
		maxRHS = finalize_RHS();

		finalize_LHS(A,b,x,1);
		finalize_LHS(A,b,x,2);
		finalize_LHS(A,b,x,3);

		assemble_RHS(b,x);

		finalize_ksp(A,b,x,1);

		for (FACE = DB.FACE; FACE; FACE = FACE->next) {
			free(FACE->LHSInIn);   FACE->LHSInIn   = NULL;
			free(FACE->LHSOutIn);  FACE->LHSOutIn  = NULL;
			free(FACE->LHSInOut);  FACE->LHSInOut  = NULL;
			free(FACE->LHSOutOut); FACE->LHSOutOut = NULL;
		}

		break;
	case 1: // diagonal VOLUME contributions
		for (VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next) {
			IndA  = VOLUME->IndA;
			NvnS  = VOLUME->NvnS;

			m = malloc(NvnS * sizeof *m); // free
			n = malloc(NvnS * sizeof *n); // free

			for (eq = 0; eq < Neq; eq++) {
				Indm = IndA + eq*NvnS;
				for (i = 0; i < NvnS; i++)
					m[i] = Indm+i;

				for (var = 0; var < Nvar; var++) {
					Indn = IndA + var*NvnS;
					for (i = 0; i < NvnS; i++)
						n[i] = Indn+i;

					vv = &(VOLUME->LHS[(eq*Nvar+var)*NvnS*NvnS]);

					MatSetValues(*A,NvnS,m,NvnS,n,vv,ADD_VALUES);
				}
			}
			free(m); free(n);
		}
		break;
	case 2: // diagonal FACE contributions
		for (FACE = DB.FACE; FACE; FACE = FACE->next) {
		for (side = 0; side < 2; side++) {
			if (side == 0) {
				VOLUME = FACE->VIn;
				LHS = FACE->LHSInIn;
			} else {
				if (FACE->Boundary)
					continue;
				VOLUME = FACE->VOut;
				LHS = FACE->LHSOutOut;
			}

			IndA = VOLUME->IndA;
			NvnS = VOLUME->NvnS;

			m = malloc(NvnS * sizeof *m); // free
			n = malloc(NvnS * sizeof *n); // free

			for (eq = 0; eq < Neq; eq++) {
				Indm = IndA + eq*NvnS;
				for (i = 0; i < NvnS; i++)
					m[i] = Indm+i;

				for (var = 0; var < Nvar; var++) {
					Indn = IndA + var*NvnS;
					for (i = 0; i < NvnS; i++)
						n[i] = Indn+i;

					vv = &LHS[(eq*Nvar+var)*NvnS*NvnS];

					MatSetValues(*A,NvnS,m,NvnS,n,vv,ADD_VALUES);
				}
			}
			free(m); free(n);
		}}
		break;
	case 3: // off-diagonal contributions
		for (FACE = DB.FACE; FACE; FACE = FACE->next) {
		for (side = 0; side < 2; side++) {
			if (FACE->Boundary)
				continue;

			if (side == 0) {
				VOLUME  = FACE->VOut;
				VOLUME2 = FACE->VIn;
				LHS = FACE->LHSInOut;
			} else {
				VOLUME  = FACE->VIn;
				VOLUME2 = FACE->VOut;
				LHS = FACE->LHSOutIn;
			}

			IndA  = VOLUME->IndA;
			IndA2 = VOLUME2->IndA;
			NvnS  = VOLUME->NvnS;
			NvnS2 = VOLUME2->NvnS;

			m = malloc(NvnS  * sizeof *m); // free
			n = malloc(NvnS2 * sizeof *n); // free

			for (eq = 0; eq < Neq; eq++) {
				Indm = IndA + eq*NvnS;
				for (i = 0; i < NvnS; i++)
					m[i] = Indm+i;

				for (var = 0; var < Nvar; var++) {
					Indn = IndA2 + var*NvnS2;
					for (i = 0; i < NvnS2; i++)
						n[i] = Indn+i;

					vv = &LHS[(eq*Nvar+var)*NvnS*NvnS2];

					MatSetValues(*A,NvnS,m,NvnS2,n,vv,ADD_VALUES);
				}
			}
			free(m); free(n);
		}}
		break;
	}
	return maxRHS;
}

void compute_dof(void)
{
	/*
	 *	Comments:
	 *		Likely put in a reordering algorithm based on conclusions of Persson (2008). (ToBeDeleted)
	 *
	 *	References:
	 *		Persson(2008)-Newton-GMRES_Preconditioning_for_Discontinuous_Galerkin_Discretizations_of_the_Navier-Stokes_
	 *		              Equations
	 */
	// Initialize DB Parameters
	unsigned int Nvar = DB.Nvar;

	// Standard datatypes
	unsigned int dof, nnz_d;

	struct S_VOLUME *VOLUME, *VIn, *VOut;
	struct S_FACE  *FACE;

	dof = 0;
	for (VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next) {
		VOLUME->IndA = dof;

		nnz_d = Nvar*(VOLUME->NvnS);
		dof += nnz_d;

		VOLUME->nnz_d = nnz_d;
		VOLUME->nnz_o = 0;
	}

	for (FACE = DB.FACE; FACE; FACE = FACE->next) {
		VIn  = FACE->VIn;
		VOut = FACE->VOut;

		if (VIn->indexg != VOut->indexg) {
			VIn->nnz_o  += VOut->nnz_d;
			VOut->nnz_o += VIn->nnz_d;
		}
	}

	DB.dof = dof;
}
