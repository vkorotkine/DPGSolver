// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#include "compute_GradW_DG.h"
#include "solver.h"

#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>

#include "Macros.h"
#include "Parameters.h"
#include "S_DB.h"
#include "S_VOLUME.h"
#include "S_FACE.h"

#include "solver_functions.h"
#include "matrix_functions.h"
#include "array_free.h"

#include "array_print.h"

/*
 *	Purpose:
 *		Compute weak gradient (auxiliary variable) contributions required for the computation of viscous fluxes.
 *
 *	Comments:
 *		For all 2nd order equations, the auxiliary variable, Qhat = QhatV + QhatF, is first computed here and then
 *		substituted into the primary equation such that a primal formulation of the DG scheme is obtained.  Thus,
 *		despite the presence of the two equations, a mixed formulation is not being used to solve the system; in that
 *		case, both What and Qhat would likely be global unknowns. Based on the mathematical analysis in Arnold(2002),
 *		restructuring the code to solve the system directly in the primal form would potentially have advantages in
 *		terms of:
 *			1) Consistency with the mathematical theory;
 *			2) Potential elimination of redundant storage/operations: (ToBeModified)
 *				- storage of Qhat, QhatV, QhatF;
 *				- computation of the same symmetric operators twice;
 *				- computation of both contributions of the lifting operator (for CDG2);
 *				- reordering operation in the penalty terms (Investigate, ToBeModified);
 *
 *	References:
 *		Arnold(2002)-Unified Analysis of Discontinuous Galerkin Methods for Elliptic Problems
 */

static void allocate_GradW         (const struct S_solver_info*const solver_info);
static void compute_GradW_VOLUME   (const char imex_type);
static void compute_GradW_FACE     (const char imex_type);
static void compute_GradW_finalize (const char imex_type);

void compute_GradW_DG (const struct S_solver_info*const solver_info)
{
	if (!DB.Viscous)
		return;

	if (solver_info->display)
		printf("G");

	allocate_GradW(solver_info);
	compute_GradW_VOLUME(solver_info->imex_type);
	compute_GradW_FACE(solver_info->imex_type);
	compute_GradW_finalize(solver_info->imex_type);
}

//static void check_NULL_and_allocate (double**const A, const size_t n)
void __attribute__ ((noinline)) check_NULL_and_allocate (double**const A, const size_t n)
{
	if (*A != NULL) {
//		printf("%f\n",*A[n]);
//		EXIT_UNSUPPORTED;
	}

	*A = calloc(n , sizeof **A);
}

static void allocate_GradW (const struct S_solver_info*const solver_info)
{
	/*
	 *	Purpose:
	 *		Allocate memory for storage of data relating to the auxiliary variable.
	 *
	 *	Comments:
	 *		If the primal formulation of the scheme is adopted in the future:
	 *			- the Qhat[L/R]/[L/R] terms would no longer need to be stored for every FACE at once;
	 *			- the QhatV terms would still need to be stored for all VOLUMEs and could only be freed after the loop
	 *			  over FACEs was completed.
	 *
	 *		Memory allocated here is freed at the end of compute_RLHS.
	 *
	 *		The memory for Qhat_What terms need not be scaled by Neq*Nvar as Qhat is linear in What.
	 */

	const unsigned int d    = DB.d,
	                   Nvar = DB.Nvar;

	// VOLUME
	for (struct S_VOLUME* VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next) {
		const unsigned int NvnS = VOLUME->NvnS;
		for (size_t dim = 0; dim < d; dim++) {
//if (VOLUME->QhatV[dim] != NULL)
//	EXIT_UNSUPPORTED;
			check_NULL_and_allocate(&VOLUME->QhatV[dim],NvnS*Nvar); // keep
			if (solver_info->imex_type == 'I')
{
//if (VOLUME->QhatV_What[dim] != NULL)
//	EXIT_UNSUPPORTED;
				check_NULL_and_allocate(&VOLUME->QhatV_What[dim],NvnS*NvnS); // keep
}
		}
	}

	// FACE
}

static void compute_GradW_VOLUME (const char imex_type)
{
	/*
	 *	Purpose:
	 *		Compute intermediate VOLUME contribution to Qhat.
	 *
	 *	Comments:
	 *		This is an intermediate contribution because the multiplication by MInv is not included.
	 *		The contribution from this function is duplicated in QhatV as the local contribution is required for the
	 *		numerical flux in the second equation.
	 *		It is currently hard-coded that GradW is of the same order as the solution.
	 *		Note, if Collocation is enable, that D_Strong includes the inverse cubature weights.
	 *
	 *	References:
	 *		Zwanenburg(2016)-Equivalence_between_the_Energy_Stable_Flux_Reconstruction_and_Discontinuous_Galerkin_Schemes
	 */

	unsigned int const d    = DB.d,
	                   Nvar = DB.Nvar;

	struct S_OPERATORS_V *OPS[2];

	struct S_VDATA *const VDATA = malloc(sizeof *VDATA); // free
	VDATA->OPS = (struct S_OPERATORS_V const *const *) OPS;

	for (size_t i = 0; i < 2; i++)
		OPS[i] = malloc(sizeof *OPS[i]); // free

	struct S_Dxyz *const DxyzInfo = malloc(sizeof *DxyzInfo); // free

	for (struct S_VOLUME *VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next) {
		init_VDATA(VDATA,VOLUME);

		unsigned int const NvnS = VDATA->OPS[0]->NvnS,
		                   NvnI = VDATA->OPS[0]->NvnI;

		DxyzInfo->Nbf = VDATA->OPS[0]->NvnS;
		DxyzInfo->Nn  = VDATA->OPS[0]->NvnI;
		DxyzInfo->D   = (double const *const *const) VDATA->OPS[0]->D_Strong;
		DxyzInfo->C   = VOLUME->C_vI;

		double const *const ChiS_vI = VDATA->OPS[0]->ChiS_vI;

		double **QhatV_What;
		if (imex_type == 'E') {
			QhatV_What = malloc(d * sizeof *QhatV_What); // free
			for (size_t dim = 0; dim < d; dim++)
				QhatV_What[dim] = malloc(NvnS*NvnS * sizeof *QhatV_What[dim]); // free
		} else {
			QhatV_What = VOLUME->QhatV_What;
		}

		for (size_t dim = 0; dim < d; dim++) {
			DxyzInfo->dim = dim;
			double *const Dxyz = compute_Dxyz_strong(DxyzInfo,d); // free

			// Note: The detJ_vI term cancels with the gradient operator (Zwanenburg(2016), eq. (B.2))
			if (DB.Collocated) { // ChiS_vI == I
				for (size_t i = 0; i < NvnS*NvnS; i++)
					QhatV_What[dim][i] = Dxyz[i];
			} else {
				mm_d(CBRM,CBT,CBNT,NvnS,NvnS,NvnI,1.0,0.0,ChiS_vI,Dxyz,QhatV_What[dim]);
			}
			free(Dxyz);

			// Compute intermediate Qhat contribution
			mm_CTN_d(NvnS,Nvar,NvnS,QhatV_What[dim],VOLUME->What,VOLUME->QhatV[dim]);

			for (size_t i = 0; i < NvnS*Nvar; i++)
				VOLUME->Qhat[dim][i] = VOLUME->QhatV[dim][i];
		}

		if (imex_type == 'E')
			array_free2_d(d,QhatV_What);
	}

	free(VDATA);
	for (size_t i = 0; i < 2; i++)
		free(OPS[i]);
	free(DxyzInfo);
}

static void compute_GradW_FACE (const char imex_type)
{
	/*
	 *	Purpose:
	 *		Compute intermediate FACE contribution to Qhat.
	 *
	 *	Comments:
	 *		L/R is used in place of In/Out (the previous convention). (ToBeDeleted)
	 *		This is an intermediate contribution because the multiplication by MInv is not included.
	 *		It is currently hard-coded that GradW is of the same order as the solution and that a central numerical flux
	 *		is used.
	 *		Note, if Collocation is enable, that I_Weak includes the inverse cubature weights.
	 */

	struct S_OPERATORS_F *OPSL[2], *OPSR[2];
	struct S_FDATA       *const FDATAL = malloc(sizeof *FDATAL), // free
	                     *const FDATAR = malloc(sizeof *FDATAR); // free
	FDATAL->OPS = (struct S_OPERATORS_F const *const *) OPSL;
	FDATAR->OPS = (struct S_OPERATORS_F const *const *) OPSR;

	struct S_NUMERICALFLUX *const NFLUXDATA = malloc(sizeof *NFLUXDATA); // free
	FDATAL->NFLUXDATA = NFLUXDATA;
	FDATAR->NFLUXDATA = NFLUXDATA;

	struct S_DATA *const DATA = malloc(sizeof *DATA); // free
	DATA->FDATAL    = FDATAL;
	DATA->FDATAR    = FDATAR;
	DATA->NFLUXDATA = NFLUXDATA;
	DATA->feature   = 'F';
	DATA->imex_type = imex_type;

	for (size_t i = 0; i < 2; i++) {
		OPSL[i]  = malloc(sizeof *OPSL[i]);  // free
		OPSR[i]  = malloc(sizeof *OPSR[i]);  // free
	}

	for (struct S_FACE *FACE = DB.FACE; FACE; FACE = FACE->next) {
		init_FDATA(FDATAL,FACE,'L',true);
		init_FDATA(FDATAR,FACE,'R',true);

		// Compute WL_fIL and WR_fIL (i.e. as seen from the (L)eft VOLUME)
		manage_solver_memory(DATA,'A','W'); // free

		coef_to_values_fI(FDATAL,'W',imex_type);
		compute_WR_fIL(FDATAR,FDATAL->W_fIL,FDATAR->W_fIL);

		// Compute numerical flux as seen from the left VOLUME
		NFLUXDATA->WL = FDATAL->W_fIL;
		NFLUXDATA->WR = FDATAR->W_fIL;
		manage_solver_memory(DATA,'A','S'); // free

		compute_numerical_solution(FDATAL,imex_type);
		add_Jacobian_scaling_FACE(FDATAL,imex_type,'Q');

		finalize_QhatF_Weak(FDATAL,FDATAR,'L','E');
		if (imex_type == 'I')
			finalize_QhatF_Weak(FDATAL,FDATAR,'L','I');

		if (!FACE->Boundary) {
			finalize_QhatF_Weak(FDATAL,FDATAR,'R','E');
			if (imex_type == 'I')
				finalize_QhatF_Weak(FDATAL,FDATAR,'R','I');
		}

		manage_solver_memory(DATA,'F','W');
		manage_solver_memory(DATA,'F','S');
	}

	for (size_t i = 0; i < 2; i++) {
		free(OPSL[i]);
		free(OPSR[i]);
	}
	free(NFLUXDATA);
	free(FDATAL);
	free(FDATAR);
	free(DATA);
}

static void finalize_Qhat(struct S_VOLUME const *const VOLUME, unsigned int const NvnS, double *const *const Qhat)
{
	unsigned int const d    = DB.d,
	                   Nvar = DB.Nvar;

	if (DB.Collocated) {
		double const *const detJV_vI = VOLUME->detJV_vI;
		for (size_t dim = 0; dim < d; dim++) {
			for (size_t var = 0; var < Nvar; var++) {
			for (size_t n = 0; n < NvnS; n++) {
				Qhat[dim][var*NvnS+n] /= detJV_vI[n];
			}}
		}
	} else {
		double *Qhat_tmp = malloc(NvnS*Nvar * sizeof *Qhat_tmp); // free
		for (size_t dim = 0; dim < d; dim++) {
			mm_CTN_d(NvnS,Nvar,NvnS,VOLUME->MInv,Qhat[dim],Qhat_tmp);
			for (size_t var = 0; var < Nvar; var++) {
			for (size_t n = 0; n < NvnS; n++) {
				Qhat[dim][var*NvnS+n] = Qhat_tmp[var*NvnS+n];
			}}
		}
		free(Qhat_tmp);
	}
}

static void finalize_Qhat_What(struct S_VOLUME const *const VOLUME, unsigned int const NRows, unsigned int const NCols,
                               bool const variable_eqvar, double *const *const Qhat_What)
{
	/*
	 *	Comments:
	 *		The variable_eqvar flag is provided as QhatV_What is constant for all equations and variables. The same is
	 *		true for Qhat_What for FACEs which are not on the boundary. In the case of boundary FACEs, while the
	 *		contribution may vary depending on the boundary condition employed.
	 */

	unsigned int const d = DB.d;

	unsigned int eqMax, varMax;

	if (variable_eqvar) {
		eqMax  = DB.Neq;
		varMax = DB.Nvar;
	} else {
		eqMax  = 1;
		varMax = 1;
	}

	if (DB.Collocated) {
		double const *const detJV_vI = VOLUME->detJV_vI;
		for (size_t dim = 0; dim < d; dim++) {
			for (size_t eq = 0; eq < eqMax; eq++) {
			for (size_t var = 0; var < varMax; var++) {
				size_t const Indeqvar = (eq*varMax+var)*NRows*NCols;
				for (size_t i = 0; i < NRows; i++) {
				for (size_t j = 0; j < NCols; j++) {
					Qhat_What[dim][Indeqvar+i*NCols+j] /= detJV_vI[i];
				}}
			}}
		}
	} else {
		double *Qhat_tmp = malloc(NRows*NCols * sizeof *Qhat_tmp); // free
		for (size_t dim = 0; dim < d; dim++) {
			for (size_t eq = 0; eq < eqMax; eq++) {
			for (size_t var = 0; var < varMax; var++) {
				size_t const Indeqvar = (eq*varMax+var)*NRows*NCols;
				mm_d(CBRM,CBNT,CBNT,NRows,NCols,NRows,1.0,0.0,VOLUME->MInv,&Qhat_What[dim][Indeqvar],Qhat_tmp);
				for (size_t i = 0; i < NRows; i++) {
				for (size_t j = 0; j < NCols; j++) {
					Qhat_What[dim][Indeqvar+i*NCols+j] = Qhat_tmp[i*NCols+j];
				}}
			}}
		}
		free(Qhat_tmp);
	}
}

static void compute_GradW_finalize (const char imex_type)
{
	/*
	 *	Purpose:
	 *		Add inverse mass matrix contribution to VOLUME->Qhat, VOLUME->QhatV.
	 *
	 *	Comments:
	 *		The contributions to Qhat_What from VOLUMEs and FACEs are not summed as they are treated in the FACE info
	 *		function.
	 *
	 *		All contributions continue to be stored individually as they must be used as such for the computation of the
	 *		viscous numerical flux. The FACE Qhat contributions are added to the VOLUME Qhat contribution to store the
	 *		entire weak gradient in VOLUME->Qhat (used for VOLUME terms).
	 *
	 *		The FACE contribution were included direcly in VL/VR->Qhat while the VOLUME contributions were stored in
	 *		QhatV. The two are now summed in VL/VR->Qhat and QhatV is retained for use in the numerical flux for the
	 *		second equation of the mixed formulation.
	 *
	 *		If Collocation is enable, only the inverse Jacobian determinant must be added here.
	 */

	unsigned int const d    = DB.d,
	                   Nvar = DB.Nvar;

	// Add FACE contributions to VOLUME->Qhat then multiply by MInv
	for (struct S_FACE *FACE = DB.FACE; FACE; FACE = FACE->next) {
		struct S_VOLUME const *const VL = FACE->VL,
		                      *const VR = FACE->VR;

		unsigned int const NvnSL = VL->NvnS,
		                   NvnSR = VR->NvnS;

		for (size_t dim = 0; dim < d; dim++) {
			for (size_t i = 0; i < NvnSL*Nvar; i++)
				VL->Qhat[dim][i] += FACE->QhatL[dim][i];

			if (!FACE->Boundary) {
				for (size_t i = 0; i < NvnSR*Nvar; i++)
					VR->Qhat[dim][i] += FACE->QhatR[dim][i];
			}
		}

		finalize_Qhat(VL,NvnSL,FACE->QhatL);
		if (FACE->Boundary) {
			if (imex_type == 'I')
				finalize_Qhat_What(VL,NvnSL,NvnSL,1,FACE->QhatL_WhatL);
		} else {
			finalize_Qhat(VR,NvnSR,FACE->QhatR);
			if (imex_type == 'I') {
				finalize_Qhat_What(VL,NvnSL,NvnSL,0,FACE->QhatL_WhatL);
				finalize_Qhat_What(VL,NvnSL,NvnSR,0,FACE->QhatL_WhatR);
				finalize_Qhat_What(VR,NvnSR,NvnSL,0,FACE->QhatR_WhatL);
				finalize_Qhat_What(VR,NvnSR,NvnSR,0,FACE->QhatR_WhatR);
			}
		}
	}

	// Multiply VOLUME Qhat terms by MInv
	for (struct S_VOLUME *VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next) {
		unsigned int const NvnS = VOLUME->NvnS;

		finalize_Qhat(VOLUME,NvnS,VOLUME->Qhat);
		finalize_Qhat(VOLUME,NvnS,VOLUME->QhatV);

		if (imex_type == 'I')
			finalize_Qhat_What(VOLUME,NvnS,NvnS,0,VOLUME->QhatV_What);
	}
}

void free_GradW_DG (const struct S_solver_info*const solver_info)
{
	/*
	 *	Purpose:
	 *		Free memory allocated for auxiliary variable terms which is still in use.
	 */

	const unsigned int d = DB.d;

	// VOLUME
	for (struct S_VOLUME* VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next) {
		for (size_t dim = 0; dim < d; dim++) {
			FREE_NULL(VOLUME->QhatV[dim]);
			if (solver_info->imex_type == 'I')
				FREE_NULL(VOLUME->QhatV_What[dim]);
		}
	}

	// FACE
}
