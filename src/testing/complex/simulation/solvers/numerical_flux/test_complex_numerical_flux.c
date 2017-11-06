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
 */

#include "test_complex_numerical_flux.h"

#include <assert.h>
#include <stdlib.h>
#include <string.h>

#include "macros.h"

#include "test_complex_numerical_flux_advection.h"

#include "complex_multiarray.h"
#include "multiarray.h"
#include "vector.h"

#include "element_solver_dg.h"
#include "face.h"
#include "face_solver.h"
#include "volume.h"
#include "volume_solver.h"

#include "const_cast.h"
#include "math_functions.h"
#include "multiarray_operator.h"
#include "numerical_flux_advection.h"
#include "operator.h"
#include "simulation.h"
#include "solve_dg.h"
#include "test_case.h"

// Static function declarations ************************************************************************************* //

/// \brief Set the function pointers in \ref Numerical_Flux_Input_c.
static void set_derived_Numerical_Flux_Input_fptrs
	(struct Numerical_Flux_Input_c* num_flux_i ///< \ref Numerical_Flux_Input_c.
	);

// Interface functions ********************************************************************************************** //

struct Numerical_Flux_Input_c* constructor_Numerical_Flux_Input_c (const struct Simulation* sim)
{
	assert(sim->test_case->solver_method_curr == 'i');

	struct Numerical_Flux_Input* num_flux_i_b = constructor_Numerical_Flux_Input(sim); // destructed.

	struct Numerical_Flux_Input_c* num_flux_i = calloc(1,sizeof *num_flux_i); // free

	memcpy(num_flux_i,num_flux_i_b,sizeof(struct Numerical_Flux_Input)); // shallow copy of the base.
	destructor_Numerical_Flux_Input(num_flux_i_b);

	set_derived_Numerical_Flux_Input_fptrs(num_flux_i);

	return num_flux_i;
}

void destructor_Numerical_Flux_Input_c (struct Numerical_Flux_Input_c* num_flux_i)
{
	free(num_flux_i);
}

struct Numerical_Flux_c* constructor_Numerical_Flux_c (const struct Numerical_Flux_Input_c* num_flux_i)
{
	assert(((num_flux_i->bv_l.s != NULL && num_flux_i->bv_l.s->layout == 'C') &&
	        (num_flux_i->bv_r.s != NULL && num_flux_i->bv_r.s->layout == 'C')) ||
	       ((num_flux_i->bv_l.g != NULL && num_flux_i->bv_l.g->layout == 'C') &&
	        (num_flux_i->bv_r.g != NULL && num_flux_i->bv_r.g->layout == 'C')));

	struct Numerical_Flux_Input* num_flux_i_b = (struct Numerical_Flux_Input*) num_flux_i;

	const int n_eq = num_flux_i_b->bv_l.n_eq;
	const ptrdiff_t n_n = ( num_flux_i->bv_l.s != NULL ? num_flux_i->bv_l.s->extents[0]
	                                                   : num_flux_i->bv_l.g->extents[0] );

	struct mutable_Numerical_Flux_c* num_flux = calloc(1,sizeof *num_flux); // free

	num_flux->nnf = constructor_zero_Multiarray_c('C',2,(ptrdiff_t[]){n_n,n_eq}); // destructed

	num_flux_i->compute_Numerical_Flux(num_flux_i,num_flux);

	return (struct Numerical_Flux_c*) num_flux;
}

void destructor_Numerical_Flux_c (struct Numerical_Flux_c* num_flux)
{
	destructor_const_Multiarray_c(num_flux->nnf);
	free(num_flux);
}

void compute_Numerical_Flux_c_1
	(const struct Numerical_Flux_Input_c* num_flux_i, struct mutable_Numerical_Flux_c* num_flux)
{
	num_flux_i->compute_Numerical_Flux_1st(num_flux_i,num_flux);
}

void compute_Numerical_Flux_c_12
	(const struct Numerical_Flux_Input_c* num_flux_i, struct mutable_Numerical_Flux_c* num_flux)
{
// make sure to use '+=' for the numerical fluxes.
	num_flux_i->compute_Numerical_Flux_1st(num_flux_i,num_flux);
	num_flux_i->compute_Numerical_Flux_2nd(num_flux_i,num_flux);
EXIT_ERROR("Ensure that all is working correctly.");
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

static void set_derived_Numerical_Flux_Input_fptrs (struct Numerical_Flux_Input_c* num_flux_i)
{
	struct Numerical_Flux_Input* num_flux_i_b = (struct Numerical_Flux_Input*) num_flux_i;

	// compute_Numerical_Flux
	if (num_flux_i_b->compute_Numerical_Flux == compute_Numerical_Flux_1)
		num_flux_i->compute_Numerical_Flux = compute_Numerical_Flux_c_1;
	else if (num_flux_i_b->compute_Numerical_Flux == compute_Numerical_Flux_12)
		num_flux_i->compute_Numerical_Flux = compute_Numerical_Flux_c_12;
	else
		EXIT_UNSUPPORTED;

	// compute_Numerical_Flux_1st
	if (num_flux_i_b->compute_Numerical_Flux_1st == compute_Numerical_Flux_advection_upwind_jacobian)
		num_flux_i->compute_Numerical_Flux_1st = compute_Numerical_Flux_c_advection_upwind;
//	else if (num_flux_i_b->compute_Numerical_Flux_1st == compute_Numerical_Flux_euler_roe_pike_jacobian)
//		num_flux_i->compute_Numerical_Flux_1st = compute_Numerical_Flux_c_euler_roe_pike;
	else if (num_flux_i_b->compute_Numerical_Flux_1st == NULL)
		num_flux_i->compute_Numerical_Flux_1st = NULL;
	else
		EXIT_UNSUPPORTED;

	// compute_Numerical_Flux_2nd
	if (num_flux_i_b->compute_Numerical_Flux_2nd == NULL)
		num_flux_i->compute_Numerical_Flux_2nd = NULL;
	else
		EXIT_UNSUPPORTED;
}
