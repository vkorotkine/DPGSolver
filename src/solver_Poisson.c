// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#include "solver_Poisson.h"

#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>

#include "petscvec.h"
#include "petscmat.h"
#include "petscksp.h"

#include "Parameters.h"
#include "Macros.h"
#include "S_DB.h"
#include "S_VOLUME.h"
#include "S_FACE.h"
#include "Test.h"

#include "solver_functions.h"
#include "update_VOLUMEs.h"
#include "implicit_GradW.h"

#include "matrix_functions.h"
#include "finalize_LHS.h"
#include "solver_implicit.h"
#include "output_to_paraview.h"
#include "test_code_output_to_paraview.h"

/*
 *	Purpose:
 *		Perform the implicit solve for the Poisson equation.
 *
 *	Comments:
 *		Many of the RHS terms computed are 0; they are included as they are used to check the linearization. Further,
 *		the computational cost is dominated by the global system solve making this additional cost negligible.
 *
 *		Modify the implementation such that the implicit_VOLUME/FACE functions can be called directly here (Analogously
 *		to solver_Advection). (ToBeModified)
 *
 *	Notation:
 *
 *	References:
 */

static void compute_Qhat(void)
{
	/*
	 *	Purpose:
	 *		Compute the weak gradients and store LHSQ for use below.
	 */

	implicit_GradW_VOLUME();

	unsigned int d = DB.d;

	for (struct S_VOLUME *VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next) {
		unsigned int const NvnS = VOLUME->NvnS;

		// Store contribution before multiplication by MInv
		for (size_t dim = 0; dim < d; dim++) {
			for (size_t i = 0; i < NvnS*NvnS; i++)
				VOLUME->LHSQ[dim][i] = -VOLUME->QhatV_What[dim][i];
			mkl_dimatcopy('R','T',NvnS,NvnS,1.0,VOLUME->LHSQ[dim],NvnS,NvnS);
		}
	}

	implicit_GradW_FACE();
	implicit_GradW_finalize();
}

static void compute_What_VOLUME(void)
{
	unsigned int const d    = DB.d,
	                   Neq  = 1,
	                   Nvar = 1;

	for (struct S_VOLUME *VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next) {
		unsigned int const NvnS = VOLUME->NvnS;

		// Compute RHS and LHS terms
		double **LHSQ = VOLUME->LHSQ;

		// RHS
		memset(VOLUME->RHS,0.0,NvnS*Neq * sizeof *(VOLUME->RHS));
		for (size_t dim = 0; dim < d; dim++)
			mm_d(CBCM,CBT,CBNT,NvnS,1,NvnS,1.0,1.0,LHSQ[dim],VOLUME->Qhat[dim],VOLUME->RHS);

		// LHS
		memset(VOLUME->LHS,0.0,NvnS*NvnS*Neq*Nvar * sizeof *(VOLUME->LHS));
		for (size_t dim = 0; dim < d; dim++)
			mm_d(CBRM,CBNT,CBNT,NvnS,NvnS,NvnS,1.0,1.0,LHSQ[dim],VOLUME->QhatV_What[dim],VOLUME->LHS);
	}
}

static void compute_What_FACE()
{
	/*
	 *	Comments:
	 *		The flux and flux Jacobians are negated below such that finalize_FACE_Inviscid_Weak can be used for the
	 *		computation of the FACE contributions.
	 */

	unsigned int d = DB.d;

	struct S_OPERATORS_F *OPSL[2], *OPSR[2];

	struct S_FDATA *FDATAL = malloc(sizeof *FDATAL), // free
	               *FDATAR = malloc(sizeof *FDATAR); // free
	FDATAL->OPS = (struct S_OPERATORS_F const *const *) OPSL;
	FDATAR->OPS = (struct S_OPERATORS_F const *const *) OPSR;

	struct S_NUMERICALFLUX *NFLUXDATA = malloc(sizeof *NFLUXDATA); // free
	FDATAL->NFLUXDATA = NFLUXDATA;
	FDATAR->NFLUXDATA = NFLUXDATA;

	struct S_DATA *const DATA = malloc(sizeof *DATA); // free
	DATA->FDATAL    = FDATAL;
	DATA->FDATAR    = FDATAR;
	DATA->NFLUXDATA = NFLUXDATA;
	DATA->feature   = 'F';
	DATA->imex_type = 'I';

	for (size_t i = 0; i < 2; i++) {
		OPSL[i] = malloc(sizeof *OPSL[i]); // free
		OPSR[i] = malloc(sizeof *OPSR[i]); // free
	}

	for (struct S_FACE *FACE = DB.FACE; FACE; FACE = FACE->next) {
		init_FDATA(FDATAL,FACE,'L');
		init_FDATA(FDATAR,FACE,'R');

		// Compute WL_fIL, WR_fIL, QpL_fIL, and QpR_fIL (i.e. as seen from the (L)eft VOLUME)
		manage_solver_memory(DATA,'A','W'); // free
		manage_solver_memory(DATA,'A','Q'); // free

		coef_to_values_fI(FDATAL,'W','I');
		coef_to_values_fI(FDATAL,'Q','I');
		compute_WR_QpR_fIL(FDATAR,FDATAL->W_fIL,FDATAR->W_fIL,
		                   (double const *const *const) FDATAL->Qp_fIL,FDATAR->Qp_fIL,'I');

		// Compute numerical flux and its Jacobian as seen from the left VOLUME
		NFLUXDATA->WL = FDATAL->W_fIL;
		NFLUXDATA->WR = FDATAR->W_fIL;
		manage_solver_memory(DATA,'A','P'); // free

		compute_numerical_flux_viscous(FDATAL,FDATAR,'I');
		manage_solver_memory(DATA,'F','W');
		manage_solver_memory(DATA,'F','Q');

		add_Jacobian_scaling_FACE(FDATAL,'E','V');
		add_Jacobian_scaling_FACE(FDATAL,'I','P');

		// Negate the flux and flux Jacobian terms
		unsigned int const IndFType = FDATAL->IndFType,
		                   NfnI     = OPSL[IndFType]->NfnI;

		for (size_t n = 0; n < NfnI; n++)
			NFLUXDATA->nFluxNum[n] *= -1.0;

		for (size_t dim = 0; dim < d; dim++) {
			for (size_t n = 0; n < NfnI; n++) {
				NFLUXDATA->dnFluxNumdQL[dim][n] *= -1.0;
				NFLUXDATA->dnFluxNumdQR[dim][n] *= -1.0;
			}
		}

		// Finalize FACE RHS and LHS terms
		unsigned int const NvnSL = OPSL[0]->NvnS,
		                   NvnSR = OPSR[0]->NvnS;

		// Add VOLUME contributions to RHS and LHS
		memset(FACE->RHSL, 0.0,NvnSL       * sizeof *FACE->RHSL);
		memset(FACE->LHSLL,0.0,NvnSL*NvnSL * sizeof *FACE->LHSLL);

		if (!FACE->Boundary) {
			memset(FACE->RHSR, 0.0,NvnSR       * sizeof *FACE->RHSR);
			memset(FACE->LHSRL,0.0,NvnSL*NvnSR * sizeof *FACE->LHSRL);
			memset(FACE->LHSLR,0.0,NvnSR*NvnSL * sizeof *FACE->LHSLR);
			memset(FACE->LHSRR,0.0,NvnSR*NvnSR * sizeof *FACE->LHSRR);
		}
		finalize_VOLUME_LHSQF_Weak(FACE);

		// Add FACE contributions to RHS and LHS

		// Interior FACE

		finalize_FACE_Viscous_Weak(FDATAL,FDATAR,NFLUXDATA->nFluxNum,NULL,'L','E','V');
		finalize_implicit_FACE_Q_Weak(FDATAL,FDATAR,'L');

		// Exterior FACE
		if (!FACE->Boundary) {
			finalize_FACE_Viscous_Weak(FDATAL,FDATAR,NFLUXDATA->nFluxNum,NULL,'R','E','V');
			finalize_implicit_FACE_Q_Weak(FDATAL,FDATAR,'R');
		}
		manage_solver_memory(DATA,'F','P');
	}
	for (size_t i = 0; i < 2; i++) {
		free(OPSL[i]);
		free(OPSR[i]);
	}

	free(FDATAL);
	free(FDATAR);
	free(NFLUXDATA);
	free(DATA);
}

void implicit_info_Poisson(void)
{
	update_VOLUME_Ops();
	update_memory_VOLUMEs();

	compute_Qhat();

	compute_What_VOLUME();
	compute_What_FACE();
}

void solver_Poisson(bool PrintEnabled)
{
	if (!strstr(DB.SolverType,"Implicit"))
		EXIT_UNSUPPORTED;

	implicit_info_Poisson();
//	if (PrintEnabled) { printf("V");  } implicit_VOLUME_info();
//	if (PrintEnabled) { printf("F");  } implicit_FACE_info();

	Mat A = NULL;
	Vec b = NULL, x = NULL;
	KSP ksp = NULL;

	solver_implicit_linear_system(&A,&b,&x,&ksp,0,PrintEnabled);
	solver_implicit_update_What(x);

	KSPDestroy(&ksp);
	finalize_ksp(&A,&b,&x,2);

	// Update Qhat based on computed solution
	implicit_GradW();

	char *const fNameOut = get_fNameOut("SolFinal_"); // free
	output_to_paraview(fNameOut);
	free(fNameOut);
}
