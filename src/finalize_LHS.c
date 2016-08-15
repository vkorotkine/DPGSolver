// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

#include "finalize_LHS.h"

#include <stdlib.h>
#include <stdio.h>

#include "petscmat.h"

#include "Macros.h"
#include "S_DB.h"
#include "S_VOLUME.h"
#include "S_FACET.h"

#include "array_print.h"

/*
 *	Purpose:
 *		Finalize LHS term by summing VOLUME and FACET contributions.
 *
 *	Comments:
 *
 *	Notation:
 *
 *	References:
 */

void initialize_KSP(Mat *A, Vec *b)
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

	*b = NULL;
	if (*b)
		printf("placeholder\n");

	free(nnz);
}

void finalize_LHS(Mat *A, Vec *b, const unsigned int assemble_type)
{
	/*
	 *	Comments:
	 *		assemble_type is a flag for which parts of the global matrix are assembled.
	 */

	// Initialize DB Parameters
	unsigned int Nvar = DB.Nvar,
	             Neq  = DB.Neq;

	// Standard datatypes
	unsigned int i, NvnS, eq, var,
	             IndA, Indm, Indn;

	struct S_VOLUME *VOLUME;

	PetscInt    *m, *n;
	PetscScalar *vv;

	compute_dof();

	if (*A == NULL)
		initialize_KSP(A,b);

	switch (assemble_type) {
	default: // 0
		printf("Error: Not yet implemented.\n"), EXIT_MSG;
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

					MatSetValues(*A,NvnS,m,NvnS,n,vv,INSERT_VALUES);
				}
			}
			free(m); free(n);
		}
		break;
	case 2: // diagonal FACET contributions
		break;
	case 3: // offdiagonal contributions
		break;
	}

	MatAssemblyBegin(*A,MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(*A,MAT_FINAL_ASSEMBLY);
//	MatView(*A,PETSC_VIEWER_STDOUT_SELF);
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
	struct S_FACET  *FACET;

	dof = 0;
	for (VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next) {
		VOLUME->IndA = dof;

		nnz_d = Nvar*(VOLUME->NvnS);
		dof += nnz_d;

		VOLUME->nnz_d = nnz_d;
		VOLUME->nnz_o = 0;
	}

	for (FACET = DB.FACET; FACET; FACET = FACET->next) {
		VIn  = FACET->VIn;
		VOut = FACET->VOut;

		if (VIn->indexg != VOut->indexg) {
			VIn->nnz_o  += VOut->nnz_d;
			VOut->nnz_o += VIn->nnz_d;
		}
	}

	DB.dof = dof;
}
