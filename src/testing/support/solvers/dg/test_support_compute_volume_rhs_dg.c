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
 *  \todo Attempt to template these functions.
 */

#include "test_support_compute_volume_rhs_dg.h"

#include <assert.h>

#include "macros.h"
#include "definitions_intrusive.h"

#include "test_support_math_functions.h"

#include "compute_volume_rlhs_dg.h"
#include "intrusive.h"
#include "simulation.h"

// Static function declarations ************************************************************************************* //

/** \brief Constructor for a \ref Flux_Ref container.
 *  \return See brief. */
/*static struct Flux_Ref* constructor_Flux_Ref
	(const struct const_Multiarray_d* m, ///< The metric terms.
	 const struct Flux* flux             ///< The physical \ref Flux.
	); */

/// \brief Compute only the rhs term.
/*static void compute_rhs
	(const struct Flux_Ref* flux_r, ///< Defined for \ref compute_rlhs_fptr.
	 struct Volume* volume,         ///< Defined for \ref compute_rlhs_fptr.
	 const struct Simulation* sim   ///< Defined for \ref compute_rlhs_fptr.
	); */

// Interface functions ********************************************************************************************** //

void compute_volume_rhs_dg_c (const struct Simulation* sim, struct Intrusive_List* volumes)
{
	assert(sim->volumes->name == IL_VOLUME_SOLVER_DG_COMPLEX);
	assert(sim->elements->name == IL_ELEMENT_SOLVER_DG);

//	struct S_Params s_params = set_s_params(sim);
//	struct Flux_Input* flux_i = constructor_Flux_Input(sim); // destructed

	for (struct Intrusive_Link* curr = volumes->first; curr; curr = curr->next) {
// Change all relevant functions to static.
UNUSED(curr);
/*		struct Volume*        vol   = (struct Volume*) curr;
		struct Solver_Volume* s_vol = (struct Solver_Volume*) curr;

		// Compute the solution, gradients and xyz coordinates at the volume cubature nodes.
		flux_i->s   = s_params.constructor_sol_vc(vol,sim);
		flux_i->g   = NULL;
		flux_i->xyz = NULL;;

		// Compute the fluxes (and optionally their Jacobians) at the volume cubature nodes.
		struct Flux* flux = constructor_Flux(flux_i);
		s_params.destructor_sol_vc(flux_i->s);

		// Compute the reference fluxes (and optionally their Jacobians) at the volume cubature nodes.
		struct Flux_Ref* flux_r = constructor_Flux_Ref(s_vol->metrics_vc,flux);
		destructor_Flux(flux);

		compute_rhs(flux_r,vol,sim);
		destructor_Flux_Ref(flux_r); */
	}
//	destructor_Flux_Input(flux_i);
EXIT_ADD_SUPPORT;
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

/** \brief Constructor for a \ref const_Multiarray_c\* of reference flux from physical flux.
 *  \return See brief. */
/*static const struct const_Multiarray_c* constructor_flux_ref
	(const struct const_Multiarray_d* m, ///< Defined for \ref constructor_Flux_Ref.
	 const struct const_Multiarray_c* f  ///< The physical flux data.
	);*/

/*static struct Flux_Ref* constructor_Flux_Ref (const struct const_Multiarray_d* m, const struct Flux* flux)
{
	assert(flux->f != NULL);
	assert(m->extents[0] == flux->f->extents[0]);

	struct Flux_Ref* flux_r = calloc(1,sizeof *flux_r); // returned

	flux_r->fr = constructor_flux_ref(m,flux->f);

	return flux_r;
}*/

/*static void compute_rhs (const struct Flux_Ref* flux_r, struct Volume* volume, const struct Simulation* sim)
{
	const struct Multiarray_Operator tw1_vs_vc = get_operator__tw1_vs_vc__rlhs_dg(volume);

	struct Complex_DG_Solver_Volume* dg_s_volume = (struct Complex_DG_Solver_Volume*) volume;
	const ptrdiff_t d = sim->d;
	for (ptrdiff_t dim = 0; dim < d; ++dim)
		mm_NNC_Operator_Multiarray_c(1.0,1.0,tw1_vs_vc.data[dim],flux_r->fr,dg_s_volume->rhs,'d',2,&dim,NULL);
}*/

// Level 1 ********************************************************************************************************** //

/*static const struct const_Multiarray_c* constructor_flux_ref
	(const struct const_Multiarray_d* m, const struct const_Multiarray_c* f)
{
	assert(f->layout == 'C');

	const int order = f->order;
	ptrdiff_t extents[order];
	for (int i = 0; i < order; ++i) {
		if (i == 0)
			extents[i] = f->extents[i];
		else if (i == 1)
			extents[order-1] = f->extents[i];
		else
			extents[i-1] = f->extents[i];
	}

	struct Multiarray_c* fr = constructor_zero_Multiarray_c('C',order,extents); // returned

	const int n_n   = extents[0];
	const int d     = extents[order-1];
	const int n_col = compute_size(order,extents)/(n_n*d);

	int ind_f = 0;
	for (int col = 0; col < n_col; ++col) {
		for (int dim0 = 0; dim0 < d; ++dim0) {
			const int ind_fr = (ind_f+dim0*n_col);
			for (int dim1 = 0; dim1 < d; ++dim1) {
				const int ind_m  = dim0*d+dim1;
				const int ind_fp = (ind_f*d)+dim1;
				z_yxpz_dcc(n_n,get_col_const_Multiarray_d(ind_m,m),
				               get_col_const_Multiarray_c(ind_fp,f),
				               get_col_Multiarray_c(ind_fr,fr));
			}
		}
		++ind_f;
	}

	return (const struct const_Multiarray_c*) fr;
}*/
