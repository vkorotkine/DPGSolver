// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#include "implicit_FACE_info_HDG.h"

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

static void compute_Inviscid_FACE_HDG (void);

void implicit_FACE_info_HDG (bool const PrintEnabled)
{
	if (PrintEnabled)
		printf("F");

	compute_Inviscid_FACE_HDG();
}

static void compute_Inviscid_FACE_HDG (void)
{
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
