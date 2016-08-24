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
printf("inKSP: %d \n",dof);

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

void finalize_Mat(Mat *A, const unsigned int finalize_type)
{
	switch (finalize_type) {
		case 1:
			MatAssemblyBegin(*A,MAT_FINAL_ASSEMBLY);
			MatAssemblyEnd(*A,MAT_FINAL_ASSEMBLY);
			break;
		default:
			printf("Error: Unsupported finalize_type.\n"), EXIT_MSG;
			break;
	}
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
	unsigned int i, NvnS, NvnS2, eq, var, side,
	             IndA, IndA2, Indm, Indn;
	double       *LHS;

	struct S_VOLUME *VOLUME, *VOLUME2;
	struct S_FACET  *FACET;

	PetscInt    *m, *n;
	PetscScalar *vv;

	compute_dof();

	if (*A == NULL)
		initialize_KSP(A,b);

	switch (assemble_type) {
	default: // 0
		finalize_LHS(A,b,1);
		finalize_LHS(A,b,2);
		finalize_LHS(A,b,3);

		finalize_Mat(A,1);

		for (FACET = DB.FACET; FACET; FACET = FACET->next) {
			free(FACET->LHSInIn);
			free(FACET->LHSOutIn);
			free(FACET->LHSInOut);
			free(FACET->LHSOutOut);
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
	case 2: // diagonal FACET contributions
		for (FACET = DB.FACET; FACET; FACET = FACET->next) {
		for (side = 0; side < 2; side++) {
			if (side == 0) {
				VOLUME = FACET->VIn;
				LHS = FACET->LHSInIn;
			} else {
				if (FACET->Boundary)
					continue;
				VOLUME = FACET->VOut;
				LHS = FACET->LHSOutOut;
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
	case 3: // offdiagonal contributions
		for (FACET = DB.FACET; FACET; FACET = FACET->next) {
		for (side = 0; side < 2; side++) {
			if (FACET->Boundary)
				continue;

			if (side == 0) {
				VOLUME  = FACET->VOut;
				VOLUME2 = FACET->VIn;
				LHS = FACET->LHSInOut;
			} else {
				VOLUME  = FACET->VIn;
				VOLUME2 = FACET->VOut;
				LHS = FACET->LHSOutIn;
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
/*
printf("%d %d %d\n",FACET->indexg,FACET->Boundary,side);
array_print_i(1,NvnS,m,'R');
array_print_i(1,NvnS2,n,'R');
*/
					MatSetValues(*A,NvnS,m,NvnS2,n,vv,ADD_VALUES);
				}
			}
			free(m); free(n);
		}}
		break;
	}
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
