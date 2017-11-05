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

#include "numerical_flux.h"

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>

#include "macros.h"

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
#include "operator.h"
#include "simulation.h"
#include "solve_dg.h"
#include "test_case.h"

// Static function declarations ************************************************************************************* //

/** \brief Add the linearized contribution from "right" volume to the linearized numerical flux terms if on a domain
 *         boundary. */
static void combine_num_flux_boundary
	(struct Numerical_Flux_Input* num_flux_i, ///< Defined for \ref compute_Numerical_Flux_fptr.
	 struct mutable_Numerical_Flux* num_flux  ///< Defined for \ref compute_Numerical_Flux_fptr.
	);

// Interface functions ********************************************************************************************** //

struct Numerical_Flux_Input* constructor_Numerical_Flux_Input (const struct Simulation* sim)
{
	struct Numerical_Flux_Input* num_flux_i = calloc(1,sizeof *num_flux_i); // returned

	const_cast_c1(&num_flux_i->bv_l.input_path,sim->input_path);

	struct Test_Case* test_case = sim->test_case;
	const_cast_b(&num_flux_i->has_1st_order,test_case->has_1st_order);
	const_cast_b(&num_flux_i->has_2nd_order,test_case->has_2nd_order);

	const_cast_i(&num_flux_i->bv_l.d,sim->d);
	const_cast_i(&num_flux_i->bv_l.n_eq,test_case->n_eq);
	const_cast_i(&num_flux_i->bv_l.n_var,test_case->n_var);

	num_flux_i->compute_Numerical_Flux = test_case->compute_Numerical_Flux;
	switch (test_case->solver_method_curr) {
	case 'e':
		num_flux_i->bv_l.compute_member = test_case->flux_comp_mem_e;
		num_flux_i->compute_Numerical_Flux_1st = test_case->compute_Numerical_Flux_e[0];
		num_flux_i->compute_Numerical_Flux_2nd = test_case->compute_Numerical_Flux_e[1];
		break;
	case 'i':
		num_flux_i->bv_l.compute_member = test_case->flux_comp_mem_i;
		num_flux_i->compute_Numerical_Flux_1st = test_case->compute_Numerical_Flux_i[0];
		num_flux_i->compute_Numerical_Flux_2nd = test_case->compute_Numerical_Flux_i[1];
		break;
	default:
		EXIT_ERROR("Unsupported: %c.\n",test_case->solver_method_curr);
		break;
	}

	return num_flux_i;
}

void destructor_Numerical_Flux_Input (struct Numerical_Flux_Input* num_flux_i)
{
	free(num_flux_i);
}

struct Numerical_Flux* constructor_Numerical_Flux (const struct Numerical_Flux_Input* num_flux_i)
{
	assert(((num_flux_i->bv_l.s != NULL && num_flux_i->bv_l.s->layout == 'C') &&
	        (num_flux_i->bv_r.s != NULL && num_flux_i->bv_r.s->layout == 'C')) ||
	       ((num_flux_i->bv_l.g != NULL && num_flux_i->bv_l.g->layout == 'C') &&
	        (num_flux_i->bv_r.g != NULL && num_flux_i->bv_r.g->layout == 'C')));

	const bool* c_m = num_flux_i->bv_l.compute_member;
	assert(c_m[0] || c_m[1] || c_m[2]);

	const int d    = num_flux_i->bv_l.d,
	          n_eq = num_flux_i->bv_l.n_eq,
	          n_vr = num_flux_i->bv_l.n_var;
	const ptrdiff_t n_n = ( num_flux_i->bv_l.s != NULL ? num_flux_i->bv_l.s->extents[0]
	                                                   : num_flux_i->bv_l.g->extents[0] );

	struct mutable_Numerical_Flux* num_flux = calloc(1,sizeof *num_flux); // returned

	num_flux->nnf = (c_m[0] ? constructor_zero_Multiarray_d('C',2,(ptrdiff_t[]){n_n,n_eq}) : NULL); // destructed
	for (int i = 0; i < 2; ++i) {
		struct m_Neigh_Info_NF* n_i = &num_flux->neigh_info[i];
		n_i->dnnf_ds = (c_m[1] ? constructor_zero_Multiarray_d('C',3,(ptrdiff_t[]){n_n,n_eq,n_vr})   : NULL); // destructed
		n_i->dnnf_dg = (c_m[2] ? constructor_zero_Multiarray_d('C',4,(ptrdiff_t[]){n_n,n_eq,n_vr,d}) : NULL); // destructed
	}

	num_flux_i->compute_Numerical_Flux(num_flux_i,num_flux);

	return (struct Numerical_Flux*) num_flux;
}

void destructor_Numerical_Flux (struct Numerical_Flux* num_flux)
{
	if (num_flux->nnf != NULL)
		destructor_const_Multiarray_d(num_flux->nnf);

	for (int i = 0; i < 2; ++i) {
		struct Neigh_Info_NF n_i = num_flux->neigh_info[i];
		if (n_i.dnnf_ds != NULL)
			destructor_const_Multiarray_d(n_i.dnnf_ds);
		if (n_i.dnnf_dg != NULL)
			destructor_const_Multiarray_d(n_i.dnnf_dg);
	}
	free(num_flux);
}

void compute_Numerical_Flux_1 (const struct Numerical_Flux_Input* num_flux_i, struct mutable_Numerical_Flux* num_flux)
{
	num_flux_i->compute_Numerical_Flux_1st(num_flux_i,num_flux);
	combine_num_flux_boundary((struct Numerical_Flux_Input*)num_flux_i,num_flux);
}

void compute_Numerical_Flux_12 (const struct Numerical_Flux_Input* num_flux_i, struct mutable_Numerical_Flux* num_flux)
{
// make sure to use '+=' for the numerical fluxes.
	num_flux_i->compute_Numerical_Flux_1st(num_flux_i,num_flux);
	num_flux_i->compute_Numerical_Flux_2nd(num_flux_i,num_flux);
EXIT_ERROR("Ensure that all is working correctly.");
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

static void combine_num_flux_boundary
	(struct Numerical_Flux_Input* num_flux_i, struct mutable_Numerical_Flux* num_flux)
{
	// Only performed for the linearized terms.
	if (!(num_flux_i->bv_l.compute_member[1] == true || num_flux_i->bv_l.compute_member[2] == true))
		return;

	// Only performed if on a boundary.
	if (num_flux_i->bv_r.ds_ds == NULL)
		return;

	struct Multiarray_d* dnnf_ds_l = num_flux->neigh_info[0].dnnf_ds;
	const struct const_Multiarray_d* dnnf_ds_r = (struct const_Multiarray_d*) num_flux->neigh_info[1].dnnf_ds,
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
			z_yxpz(n_n,get_col_const_Multiarray_d(ind_dnnf_ds_r,dnnf_ds_r),
			           get_col_const_Multiarray_d(ind_ds_ds,ds_ds),
			           get_col_Multiarray_d(ind_dnnf_ds_l,dnnf_ds_l));
		}
	}}

	destructor_const_Multiarray_d(num_flux_i->bv_r.ds_ds);
	num_flux_i->bv_r.ds_ds = NULL;

	destructor_Multiarray_d(num_flux->neigh_info[1].dnnf_ds);
	num_flux->neigh_info[1].dnnf_ds = NULL;

	assert(num_flux_i->bv_l.compute_member[2] == false); // Add support for this in future.
}
