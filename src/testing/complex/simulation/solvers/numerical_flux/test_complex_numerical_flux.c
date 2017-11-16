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
#include "test_support_math_functions.h"

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

/// \brief `complex` version of \ref combine_num_flux_boundary.
static void combine_num_flux_boundary_c
	(struct Numerical_Flux_Input_c* num_flux_i, ///< See brief.
	 struct mutable_Numerical_Flux_c* num_flux  ///< See brief.
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

	const int d    = num_flux_i_b->bv_l.d,
	          n_eq = num_flux_i_b->bv_l.n_eq,
	          n_vr = num_flux_i_b->bv_l.n_var;
	const ptrdiff_t n_n = ( num_flux_i->bv_l.s != NULL ? num_flux_i->bv_l.s->extents[0]
	                                                   : num_flux_i->bv_l.g->extents[0] );

	struct mutable_Numerical_Flux_c* num_flux = calloc(1,sizeof *num_flux); // free

	num_flux->nnf = constructor_zero_Multiarray_c('C',2,(ptrdiff_t[]){n_n,n_eq}); // destructed
	for (int i = 0; i < 2; ++i) {
		struct m_Neigh_Info_NF_c* n_i = &num_flux->neigh_info[i];
		n_i->dnnf_ds = constructor_zero_Multiarray_c('C',3,(ptrdiff_t[]){n_n,n_eq,n_vr});   // destructed
		n_i->dnnf_dg = constructor_zero_Multiarray_c('C',4,(ptrdiff_t[]){n_n,n_eq,n_vr,d}); // destructed
	}

	num_flux_i->compute_Numerical_Flux(num_flux_i,num_flux);

	return (struct Numerical_Flux_c*) num_flux;
}

void destructor_Numerical_Flux_c (struct Numerical_Flux_c* num_flux)
{
	destructor_const_Multiarray_c(num_flux->nnf);
	for (int i = 0; i < 2; ++i) {
		struct Neigh_Info_NF_c n_i = num_flux->neigh_info[i];
		destructor_const_Multiarray_c(n_i.dnnf_ds);
		destructor_const_Multiarray_c(n_i.dnnf_dg);
	}
	free(num_flux);
}

void compute_Numerical_Flux_c_1
	(const struct Numerical_Flux_Input_c* num_flux_i, struct mutable_Numerical_Flux_c* num_flux)
{
	num_flux_i->compute_Numerical_Flux_1st(num_flux_i,num_flux);
	combine_num_flux_boundary_c((struct Numerical_Flux_Input_c*)num_flux_i,num_flux);
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
EXIT_ERROR("Add a flag for which function to use depending on sim->method.");
	if (num_flux_i_b->compute_Numerical_Flux_1st == compute_Numerical_Flux_advection_upwind_jacobian)
		num_flux_i->compute_Numerical_Flux_1st = compute_Numerical_Flux_c_advection_upwind_jacobian;
//	else if (num_flux_i_b->compute_Numerical_Flux_1st == compute_Numerical_Flux_euler_roe_pike_jacobian)
//		num_flux_i->compute_Numerical_Flux_1st = compute_Numerical_Flux_c_euler_roe_pike_jacobian;
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

static void combine_num_flux_boundary_c
	(struct Numerical_Flux_Input_c* num_flux_i, struct mutable_Numerical_Flux_c* num_flux)
{
	struct Numerical_Flux_Input* num_flux_i_b = (struct Numerical_Flux_Input*) num_flux_i;

	// Only performed for the linearized terms.
	if (!(num_flux_i_b->bv_l.compute_member[1] == true || num_flux_i_b->bv_l.compute_member[2] == true))
		return;

	// Only performed if on a boundary.
	if (num_flux_i_b->bv_r.ds_ds == NULL)
		return;

	struct Multiarray_c* dnnf_ds_l = num_flux->neigh_info[0].dnnf_ds;
	const struct const_Multiarray_c* dnnf_ds_r = (struct const_Multiarray_c*) num_flux->neigh_info[1].dnnf_ds,
	                               * ds_ds     = num_flux_i->bv_r.ds_ds;

	const int n_n  = dnnf_ds_l->extents[0],
	          n_eq = dnnf_ds_l->extents[1],
	          n_vr = dnnf_ds_l->extents[2];

	for (int eq = 0; eq < n_eq; ++eq) {
	for (int vr_l = 0; vr_l < n_vr; ++vr_l) {
		const int ind_dnnf_ds_l = eq+vr_l*n_eq;
		for (int vr_r = 0; vr_r < n_vr; ++vr_r) {
			const int ind_dnnf_ds_r = eq+vr_r*n_eq;
			const int ind_ds_ds     = vr_l+vr_r*n_vr;
			z_yxpz_ccc(n_n,get_col_const_Multiarray_c(ind_dnnf_ds_r,dnnf_ds_r),
			           get_col_const_Multiarray_c(ind_ds_ds,ds_ds),
			           get_col_Multiarray_c(ind_dnnf_ds_l,dnnf_ds_l));
		}
	}}

	destructor_const_Multiarray_c(num_flux_i->bv_r.ds_ds);
	num_flux_i->bv_r.ds_ds = NULL;

	destructor_Multiarray_c(num_flux->neigh_info[1].dnnf_ds);
	num_flux->neigh_info[1].dnnf_ds = NULL;

	assert(num_flux_i_b->bv_l.compute_member[2] == false); // Add support for this in future.
}
