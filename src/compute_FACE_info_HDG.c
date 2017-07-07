// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#include "compute_FACE_info_HDG.h"
#include "compute_RLHS.h"

#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>

#include "S_DB.h"
#include "S_FACE.h"

/*
#include <stdio.h>
#include <string.h>

#include "Parameters.h"
#include "Macros.h"
#include "S_DB.h"

#include "solver_functions.h"
#include "array_free.h"
#include "support.h"
*/

/*
 *	Purpose:
 *		Evaluate the FACE contributions to the RHS and LHS terms for the HDG scheme.
 *
 *	Comments:
 *		As compared to DG, this scheme is different in that a numerical trace unknown is now present. As a result, the
 *		contribution of the solution can be statically condensed out of the global system (See Nguyen(2015) for the
 *		idea):
 *
 *			Given the system of equations with solution unknowns, u, and trace unknowns, uT:
 *
 *				[ A B ] [u ] = [f]
 *				[ C D ] [uT] = [g]
 *
 *			the contributions of u are statically condensed out as:
 *
 *				u = inv(A)*(-B*uT+f)
 *
 *			leading to
 *
 *				[K] [uT] = [r]
 *
 *				K = -C*inv(A)*B + D
 *				r = -C*inv(A)*f + g
 *
 *		The A contribution to the system is stored in the appropriate VOLUME->LHS term. Note that A has a contribution
 *		from implicit_VOLUME_info_HDG as well.
 *
 *	References:
 *		Nguyen(2015)-A Class of Embedded Discontinuous Galerkin Methods for Computation Fluid Dynamics
 */

static void compute_Inviscid_FACE_HDG (const char imex_type);

void compute_FACE_info_HDG (const struct S_RLHS_info*const RLHS_info)
{
	if (RLHS_info->PrintEnabled)
		printf("F");

	compute_Inviscid_FACE_HDG(RLHS_info->imex_type);
}

static void compute_Inviscid_FACE_HDG (const char imex_type)
{
if (0)
printf("%d\n",imex_type);
/*
	Step 1 (Similar to DG)
	Loop over FACEs:   Complete A (store in VOLUME->LHS)
	Loop over VOLUMEs: Compute inv(A)
	Loop over FACEs:   Compute B, include premultiplication by inv(A)
	Finalize RHS (solution): Include premultiplication by inv(A)

	Step 2
	Create FACE 'connectivity' array:
		1) Initialize arrays:
			- GFToV:  size_t            of size NGFx2 (Indices of all neighbouring FACEs)
			- GF_ptr: struct S_FACE *   of size NGFx1 (Pointers to all FACEs)
		2) Sort GFToV by row and then by column

		3) Initialize empty array:
			- GFToGF: struct S_FACE *   of size NGF x 2*NFMAX*NSUBFMAX (Pointers to all related FACEs)
		4) Loop over FACEs and
			- Find indices of adjacent VOLUMEs
			- Find all FACEs connected to each VOLUME
			- Store in array

	Step 3
	Initialize Petsc Mat for the statically condensed global system:
		- dof computed based on FACE connectivity and their orders.

	Loop over FACEs:
		- Compute -C*inv(A)*B using information from GFToGF to find all non-zero components.
		- Also store diagonal FACE contributions.


	Testing:
	1) test_integration_linearization:
		- No static condensation:
			1.1) Analytical Jacobian vs. complex step for: A (Only VOLUME, and VOLUME + FACE)
			1.2)                                           A, B
			1.3)                                           A, B, C
			1.4)                                           A, B, C, D
		- With static condensation:
			1.5) Compare Petsc sparse MatMatMult -C*inv(A)*B+D with analytically computed condensed Jacobian
			1.6) Analytical Jacobian vs. complex step for condensed Jacobian.
			1.7)   "                                                    "     using hp adapted mesh.
	2) test_integration_conv_order:
		Advection default/Peterson, Euler
*/

	// Finalize contributions to A
	for (struct S_FACE *FACE = DB.FACE; FACE; FACE = FACE->next) {
		// Interpolate to FACE cubature nodes
		// Compute the stabilization matrix (store for later use)
		// Multiply by FACE Jacobian term
		// Add to VOLUME->LHS

		/*
		 *	For the upwind stabilization (linear advection), we have:
		 *		n (dot) fNum = Average component      + Stabilization component
		 *		             = 0.5*(n (dot) b)*(u+uT) + 0.5*|n (dot) b|*(u-uT)
		 *		             = n (dot) b * uT + S *(u - uT)
		 *
		 *	where
		 *
		 *		S = 0.5*(n (dot) b + |n (dot) b|)
		 *
		 *	such that we obtain a pure upwind numerical flux. Contributions from the average and stabilization
		 *	components are not combined as this is no longer possible when considering nonlinear fluxes.
		 */
	}
}
