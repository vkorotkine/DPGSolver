// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#include "implicit_VOLUME_info_HDG.h"

#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>

#include "Macros.h"
#include "S_DB.h"
#include "S_VOLUME.h"
/*
#include <string.h>

#include "Parameters.h"
#include "Macros.h"

#include "solver_functions.h"
#include "fluxes_structs.h"
#include "fluxes_inviscid.h"
#include "fluxes_viscous.h"
#include "jacobian_fluxes_inviscid.h"
#include "jacobian_fluxes_viscous.h"
#include "array_free.h"
#include "array_print.h"
*/

#include "solver_functions_HDG.h"

#include "array_print.h"

/*
 *	Purpose:
 *		Evaluate the VOLUME contributions to the RHS and LHS terms for the HDG scheme.
 *
 *	Comments:
 *		In preparation for moving to storing all data in vector/matrix struct format, necessary arrays are moved to the
 *		appropriate format here. To avoid memory leaks and remain consistent with the current code, the struct templates
 *		are then freed while retaining the array data. This can be removed as the conversion to the new format is
 *		implemented (ToBeDeleted).
 *
 *	Notation:
 *
 *	References:
 */

static void compute_Inviscid_VOLUME_HDG (void);

void implicit_VOLUME_info_HDG (bool const PrintEnabled)
{
	if (PrintEnabled)
		printf("V");

	compute_Inviscid_VOLUME_HDG();
}

static void compute_Inviscid_VOLUME_HDG (void)
{

	for (struct S_VOLUME *VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next) {
		// convert c arrays to matrix structs
		convert_to_multiarray_V(VOLUME,'A');

		struct S_VDATA VDATA;
		VDATA.OPS    = init_mat_ops_VOLUME(VOLUME);
		VDATA.VOLUME = VOLUME;

		struct S_FLUX_MA FLUXDATA;

		// Obtain W_vI
		coef_to_values_vI_MA(&VDATA,'W','A'); // free

		// Compute Flux and its Jacobian in reference space
		compute_flux_inviscid_MA(&VDATA,&FLUXDATA,'I','A'); // free
		coef_to_values_vI_MA(&VDATA,'W','F');

		// Convert to reference space
		compute_flux_ref_MA(VOLUME->C_vI_MA,FLUXDATA.F,&FLUXDATA.Fr,'A');       // free
		compute_flux_ref_MA(VOLUME->C_vI_MA,FLUXDATA.dFdW,&FLUXDATA.dFrdW,'A'); // free
		compute_flux_inviscid_MA(&VDATA,&FLUXDATA,'I','F');

		// Compute RHS and LHS terms

		// RHS
		set_to_zero_multiarray(VOLUME->RHS_MA);
		finalize_VOLUME_Inviscid_Weak_MA(FLUXDATA.Fr,VOLUME->RHS_MA,'E',&VDATA);

multiarray_print(VOLUME->RHS_MA);
EXIT_UNSUPPORTED;
		// LHS

//		compute_flux_ref_M(VOLUME->C_vI_M,FLUXDATA.F,&FLUXDATA.Fr,'F');
//		compute_flux_ref_M(VOLUME->C_vI_M,FLUXDATA.dFdW,&FLUXDATA.dFrdW,'F');

		free((void *) VDATA.OPS);

		convert_to_multiarray_V(VOLUME,'F');
	}
}
