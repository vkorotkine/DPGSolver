// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#include "compute_FACE_RLHS_DG.h"
#include "solver.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "Parameters.h"
#include "Macros.h"
#include "S_DB.h"
#include "S_FACE.h"

#include "solver_functions.h"
#include "compute_GradW_DG.h"
#include "array_free.h"
#include "support.h"


/*
 *	Purpose:
 *		Evaluate the FACE contributions to the RHS and LHS terms.
 *
 *	Comments:
 *
 *		When adaptation is enabled, the method used here to lift the FACE information to the VOLUME remains quite
 *		similar to that of the conforming case. The normal numerical flux is computed on each of the FACEs of the mesh
 *		and is used directly to form RHSL/RHSR terms to be added to the appropriate VOLUME. This is in contrast to the
 *		traditional mortar element method (based on my current understanding), which first uses an L2 projection of the
 *		normal numerical flux of all non-conforming FACEs to standard VOLUME FACEs and then computes RHSL/RHSR exactly
 *		as if the mesh were conforming. As no L2 projection is required for the current approach, the same functions
 *		can in fact be called even for the case of non-conforming discretizations allowing for:
 *			1) Reduced cost because the interpolation from FACE solution to FACE cubature nodes is not required;
 *			2) Reduced aliasing because the normal numerical flux can be computed at the cubature nodes.
 *
 *		Based on this discussion, it is unclear why this approach is not adopted instead of the traditional mortar
 *		method (Kopriva(1996)) as this alternative seems to satisfy both the conservation and outflow condition
 *		requirements which motivated the use of the mortar element method.
 *
 *		Vectorization is more involved for FACE terms as there are many more possible combinations than for VOLUMEs.
 *		Given that the VOLUME vectorization seems not to have a significant impact on performance, the FACE
 *		vectorization may not be pursued. (ToBeModified)
 *
 *	Notation:
 *		Qp : (p)artially corrected weak solution gradients (Q).
 *
 *	References:
 *		Kopriva(1996)-A Conservative Staggered-Grid Chebyshev Multidomain Method for Compressible Flows II. A Semi-Structured Method
 */

static void set_memory_to_zero_FACEs  (void);
static void compute_Inviscid_FACE_EFE (const struct S_solver_info*const solver_info);
static void compute_Viscous_FACE_EFE  (const struct S_solver_info*const solver_info);

void compute_FACE_RLHS_DG (const struct S_solver_info*const solver_info)
{
	/*
	 *	Comments:
	 *		Should only be called through compute_RLHS (solver.c).
	 */

	if (!solver_info->compute_F)
		return;

	if (solver_info->display)
		printf("F");

	set_memory_to_zero_FACEs();
	compute_Inviscid_FACE_EFE(solver_info);
	compute_Viscous_FACE_EFE(solver_info);
}

static void set_memory_to_zero_FACEs (void)
{
	unsigned int const Neq  = DB.Neq;

	for (struct S_FACE *FACE = DB.FACE; FACE; FACE = FACE->next) {
		unsigned int const NvnSL = FACE->VL->NvnS;
		set_to_zero_d(NvnSL*Neq,FACE->RHSL);

		if (!FACE->Boundary) {
			unsigned int const NvnSR = FACE->VR->NvnS;
			set_to_zero_d(NvnSR*Neq,FACE->RHSR);
		}
	}
}

static void compute_Inviscid_FACE_EFE (const struct S_solver_info*const solver_info)
{
	if (!DB.Inviscid)
		return;

	const char imex_type = solver_info->imex_type;

	struct S_OPERATORS_F *OPSL[2], *OPSR[2];

	struct S_FDATA *FDATAL = malloc(sizeof *FDATAL), // free
	               *FDATAR = malloc(sizeof *FDATAR); // free
	FDATAL->OPS = (struct S_OPERATORS_F const *const *) OPSL;
	FDATAR->OPS = (struct S_OPERATORS_F const *const *) OPSR;

	struct S_NUMERICALFLUX *NFLUXDATA = malloc(sizeof *NFLUXDATA); // free
	FDATAL->NFLUXDATA = NFLUXDATA;
	FDATAR->NFLUXDATA = NFLUXDATA;

	for (size_t i = 0; i < 2; i++) {
		OPSL[i] = malloc(sizeof *OPSL[i]); // free
		OPSR[i] = malloc(sizeof *OPSR[i]); // free
	}

	struct S_DATA *const DATA = malloc(sizeof *DATA); // free
	DATA->FDATAL    = FDATAL;
	DATA->FDATAR    = FDATAR;
	DATA->NFLUXDATA = NFLUXDATA;
	DATA->feature   = 'F';
	DATA->imex_type = imex_type;

	if (strstr(DB.Form,"Weak")) {
		for (struct S_FACE *FACE = DB.FACE; FACE; FACE = FACE->next) {
			init_FDATA(FDATAL,FACE,'L',true);
			init_FDATA(FDATAR,FACE,'R',true);

			// Compute WL_fIL and WR_fIL (i.e. as seen from the (L)eft VOLUME)
			manage_solver_memory(DATA,'A','W'); // free

			coef_to_values_fI(FDATAL,'W',imex_type);
			compute_WR_fIL(FDATAR,FDATAL->W_fIL,FDATAR->W_fIL);


			// Compute numerical flux and its Jacobian as seen from the left VOLUME
			NFLUXDATA->WL = FDATAL->W_fIL;
			NFLUXDATA->WR = FDATAR->W_fIL;
			manage_solver_memory(DATA,'A','I'); // free

			compute_numerical_flux(FDATAL,imex_type);
			add_Jacobian_scaling_FACE(FDATAL,imex_type,'W');

			manage_solver_memory(DATA,'F','W');


			// Compute FACE RHS and LHS terms

			finalize_FACE_Inviscid_Weak(FDATAL,FDATAR,NFLUXDATA->nFluxNum,NULL,'L','E','W');
			if (imex_type == 'I') {
				manage_solver_memory(DATA,'A','L'); // free
				finalize_FACE_Inviscid_Weak(FDATAL,FDATAR,NFLUXDATA->dnFluxNumdWL,NFLUXDATA->dnFluxNumdWR,'L','I','W');
			}

			if (!FACE->Boundary) {
				finalize_FACE_Inviscid_Weak(FDATAL,FDATAR,NFLUXDATA->nFluxNum,NULL,'R','E','W');
				if (imex_type == 'I') {
					finalize_FACE_Inviscid_Weak(FDATAL,FDATAR,NFLUXDATA->dnFluxNumdWL,NFLUXDATA->dnFluxNumdWR,'R','I','W');
				}
			}

			manage_solver_memory(DATA,'F','I');

			if (imex_type == 'I') {
				struct S_LHS_info LHS_info_LL = constructor_LHS_info(FDATAL->LHSL,FACE->VL,FACE->VL,ADD_VALUES,true);
				fill_PetscMat(solver_info,&LHS_info_LL);

				if (!FACE->Boundary) {
					struct S_LHS_info LHS_info_LR = constructor_LHS_info(FDATAL->LHSR,FACE->VR,FACE->VL,ADD_VALUES,true);
					fill_PetscMat(solver_info,&LHS_info_LR);
					struct S_LHS_info LHS_info_RL = constructor_LHS_info(FDATAR->LHSL,FACE->VL,FACE->VR,ADD_VALUES,true);
					fill_PetscMat(solver_info,&LHS_info_RL);
					struct S_LHS_info LHS_info_RR = constructor_LHS_info(FDATAR->LHSR,FACE->VR,FACE->VR,ADD_VALUES,true);
					fill_PetscMat(solver_info,&LHS_info_RR);
				}
				manage_solver_memory(DATA,'F','L');
			}
		}
	} else if (strstr(DB.Form,"Strong")) {
		EXIT_UNSUPPORTED;
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

static void compute_Viscous_FACE_EFE (const struct S_solver_info*const solver_info)
{
	if (!DB.Viscous)
		return;

	const char imex_type = solver_info->imex_type;

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
	DATA->imex_type = imex_type;

	for (size_t i = 0; i < 2; i++) {
		OPSL[i] = malloc(sizeof *OPSL[i]); // free
		OPSR[i] = malloc(sizeof *OPSR[i]); // free
	}

	if (strstr(DB.Form,"Weak")) {
		for (struct S_FACE *FACE = DB.FACE; FACE; FACE = FACE->next) {
			init_FDATA(FDATAL,FACE,'L',true);
			init_FDATA(FDATAR,FACE,'R',true);

			// Compute WL_fIL, WR_fIL, QpL_fIL, and QpR_fIL (i.e. as seen from the (L)eft VOLUME)
			manage_solver_memory(DATA,'A','W'); // free
			manage_solver_memory(DATA,'A','Q'); // free

			coef_to_values_fI(FDATAL,'W',imex_type);
			coef_to_values_fI(FDATAL,'Q',imex_type);
			compute_WR_QpR_fIL(FDATAR,FDATAL->W_fIL,FDATAR->W_fIL,
			                   (double const *const *const) FDATAL->Qp_fIL,FDATAR->Qp_fIL,imex_type);


			// Compute numerical flux and its Jacobian as seen from the left VOLUME
			NFLUXDATA->WL = FDATAL->W_fIL;
			NFLUXDATA->WR = FDATAR->W_fIL;
			manage_solver_memory(DATA,'A','V'); // free

			compute_numerical_flux_viscous(FDATAL,FDATAR,imex_type);
			manage_solver_memory(DATA,'F','W');
			manage_solver_memory(DATA,'F','Q');

			add_Jacobian_scaling_FACE(FDATAL,imex_type,'V');
			if (imex_type == 'I')
				add_Jacobian_scaling_FACE(FDATAL,'I','P');

			// Compute FACE RHS and LHS terms
			finalize_FACE_Viscous_Weak(FDATAL,FDATAR,NFLUXDATA->nFluxNum,NULL,'L','E','V');
			if (imex_type == 'I') {
				manage_solver_memory(DATA,'A','L'); // free
				finalize_FACE_Viscous_Weak(FDATAL,FDATAR,NFLUXDATA->dnFluxNumdWL,NFLUXDATA->dnFluxNumdWR,'L','I','V');
				finalize_implicit_FACE_Q_Weak(FDATAL,FDATAR,'L');
			}

			if (!FACE->Boundary) {
				finalize_FACE_Viscous_Weak(FDATAL,FDATAR,NFLUXDATA->nFluxNum,NULL,'R','E','V');
				if (imex_type == 'I') {
					finalize_FACE_Viscous_Weak(FDATAL,FDATAR,NFLUXDATA->dnFluxNumdWL,NFLUXDATA->dnFluxNumdWR,'R','I','V');
					finalize_implicit_FACE_Q_Weak(FDATAL,FDATAR,'R');
				}
			}
			manage_solver_memory(DATA,'F','V');

			if (imex_type == 'I') {
				struct S_LHS_info LHS_info_LL = constructor_LHS_info(FDATAL->LHSL,FACE->VL,FACE->VL,ADD_VALUES,true);
				fill_PetscMat(solver_info,&LHS_info_LL);

				if (!FACE->Boundary) {
					struct S_LHS_info LHS_info_LR = constructor_LHS_info(FDATAL->LHSR,FACE->VR,FACE->VL,ADD_VALUES,true);
					fill_PetscMat(solver_info,&LHS_info_LR);
					struct S_LHS_info LHS_info_RL = constructor_LHS_info(FDATAR->LHSL,FACE->VL,FACE->VR,ADD_VALUES,true);
					fill_PetscMat(solver_info,&LHS_info_RL);
					struct S_LHS_info LHS_info_RR = constructor_LHS_info(FDATAR->LHSR,FACE->VR,FACE->VR,ADD_VALUES,true);
					fill_PetscMat(solver_info,&LHS_info_RR);
				}
			}

			if (imex_type == 'I')
				manage_solver_memory(DATA,'F','L');
		}
	} else if (strstr(DB.Form,"Strong")) {
		EXIT_UNSUPPORTED;
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
