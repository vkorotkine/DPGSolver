/* {{{
This file is part of DPGSolver.

DPGSolver is free software: you can redistribute it and/or modify it under the terms of the GNU
General Public License as published by the Free Software Foundation, either version 3 of the
License, or any later version.

DPGSolver is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
Public License for more details.

You should have received a copy of the GNU General Public License along with DPGSolver.  If not, see
<http://www.gnu.org/licenses/>.
}}} */
/** \file
 *  \brief Provides the templated linear advection numerical flux functions.
 */

#include <assert.h>
#include <stddef.h>
#include <math.h>

#include "macros.h"
#include "definitions_bc.h"
#include "definitions_core.h"
#include "definitions_tol.h"

#include "def_templates_multiarray.h"

#include "def_templates_boundary.h"
#include "def_templates_flux.h"
#include "def_templates_math_functions.h"
#include "def_templates_numerical_flux.h"

// Static function declarations ************************************************************************************* //

#define NEQ  NEQ_BURGERS  ///< Number of equations.
#define NVAR NVAR_BURGERS ///< Number of variables.

// Interface functions ********************************************************************************************** //

void compute_Numerical_Flux_T_burgers_inviscid_lax_friedrichs
	(const struct Numerical_Flux_Input_T* num_flux_i, struct mutable_Numerical_Flux_T* num_flux)
{
	assert(DIM == 1); // Hard-coded for 1D below. Add support.

	const ptrdiff_t n_total = num_flux_i->bv_l.s->extents[0];

	const struct const_Multiarray_T*const normals = num_flux_i->bv_l.normals;

	const struct const_Multiarray_T*const u_l = num_flux_i->bv_l.s;
	const struct const_Multiarray_T*const u_r = num_flux_i->bv_r.s;

	Type*const nnf = num_flux->nnf->data;

	struct Flux_Input_T*const flux_i = constructor_Flux_Input_T(num_flux_i->sim); // destructed

	flux_i->s = u_l;
	struct Flux_T*const flux_l = constructor_Flux_T(flux_i); // destructed

	flux_i->s = u_r;
	struct Flux_T*const flux_r = constructor_Flux_T(flux_i); // destructed
	destructor_Flux_Input_T(flux_i);

	for (int n = 0; n < n_total; n++) {
		const Type ul_n = u_l->data[n];
		const Type ur_n = u_r->data[n];

		const Type fl_n = flux_l->f->data[n];
		const Type fr_n = flux_r->f->data[n];

		const Type*const normals_n = get_row_const_Multiarray_T(n,normals);
		const Type u_max_abs = max_abs_real_T(ul_n,ur_n);

		nnf[n] = 0.5*normals_n[0]*((fl_n+fr_n) + u_max_abs*(ul_n-ur_n));
	}
	destructor_Flux_T(flux_l);
	destructor_Flux_T(flux_r);
}

void compute_Numerical_Flux_T_burgers_inviscid_lax_friedrichs_jacobian
	(const struct Numerical_Flux_Input_T* num_flux_i, struct mutable_Numerical_Flux_T* num_flux
	)
{
	EXIT_ADD_SUPPORT; UNUSED(num_flux_i); UNUSED(num_flux);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

#include "undef_templates_multiarray.h"

#include "undef_templates_boundary.h"
#include "undef_templates_flux.h"
#include "undef_templates_math_functions.h"
#include "undef_templates_numerical_flux.h"
#include "undef_templates_solution_advection.h"
