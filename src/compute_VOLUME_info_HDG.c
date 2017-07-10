// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#include "compute_VOLUME_info_HDG.h"

#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>

#include "Macros.h"
#include "S_DB.h"
#include "S_VOLUME.h"

#include "solver.h"
#include "solver_functions_HDG.h"

/*
 *	Purpose:
 *		Evaluate the VOLUME contributions to the RHS and (optionally) the LHS terms for the HDG scheme.
 *
 *	Comments:
 *		The VOLUME contribution is identical to that of the DG scheme.
 *
 *		When imex_type == 'I', LHS (linearization) terms are computed as well.
 *
 *		In preparation for moving to storing all data in vector/matrix struct format, necessary arrays are moved to the
 *		appropriate format here. To avoid memory leaks and remain consistent with the current code, the struct templates
 *		are then freed while retaining the array data. This can be removed as the conversion to the new format is
 *		implemented (ToBeDeleted).
 *
 *	Notation:
 *		imex_type: Flag indicating whether the solver is being run in (im)plicit or (ex)plicit mode.
 */

static void compute_Inviscid_VOLUME_HDG (const char imex_type);

void compute_VOLUME_info_HDG (const struct S_solver_info*const solver_info)
{
	if (solver_info->display)
		printf("V");

	compute_Inviscid_VOLUME_HDG(solver_info->imex_type);
}

static void compute_Inviscid_VOLUME_HDG (const char imex_type)
{
	for (struct S_VOLUME *VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next) {
		convert_to_multiarray_V(VOLUME,'A'); // free ('F')

		struct S_VDATA VDATA;
		set_VDATA(&VDATA,VOLUME,'A'); // free ('F')

		struct S_FLUX_MA FLUXDATA;
		set_FLUXDATA(&FLUXDATA);

		// Obtain W_vI
		coef_to_values_vI_MA(&VDATA,'W','A'); // free ('F')

		// Compute Flux and its Jacobian in reference space
		compute_flux_inviscid_MA(&VDATA,&FLUXDATA,imex_type,'A'); // free ('F')
		coef_to_values_vI_MA(&VDATA,'W','F');

		// Convert to reference space
		compute_flux_ref_MA(VOLUME->C_vI_MA,FLUXDATA.F,&FLUXDATA.Fr,'A'); // free ('F')
		if (imex_type == 'I')
			compute_flux_ref_MA(VOLUME->C_vI_MA,FLUXDATA.dFdW,&FLUXDATA.dFrdW,'A'); // free ('F')
		compute_flux_inviscid_MA(&VDATA,&FLUXDATA,imex_type,'F');

		// Compute RHS and LHS terms

		// RHS
		set_to_zero_multiarray(VOLUME->RHS_MA);
		finalize_VOLUME_Inviscid_Weak_MA(FLUXDATA.Fr,VOLUME->RHS_MA,'E',&VDATA);
		compute_flux_ref_MA(VOLUME->C_vI_MA,FLUXDATA.F,&FLUXDATA.Fr,'F');

		// LHS
		if (imex_type == 'I') {
			set_to_zero_multiarray(VOLUME->LHS_L_MA);
			finalize_VOLUME_Inviscid_Weak_MA(FLUXDATA.dFrdW,VOLUME->LHS_L_MA,imex_type,&VDATA);
			compute_flux_ref_MA(VOLUME->C_vI_MA,FLUXDATA.dFdW,&FLUXDATA.dFrdW,'F');
		}

		set_VDATA(&VDATA,VOLUME,'F');

		convert_to_multiarray_V(VOLUME,'F');
	}
}
