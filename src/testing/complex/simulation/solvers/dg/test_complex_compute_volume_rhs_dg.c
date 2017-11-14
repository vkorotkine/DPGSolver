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

#include "test_complex_compute_volume_rhs_dg.h"

#include <assert.h>

#include "macros.h"
#include "definitions_intrusive.h"

#include "test_complex_flux.h"
#include "test_complex_operators.h"
#include "test_support_math_functions.h"

#include "volume_solver_dg_complex.h"

#include "complex_multiarray.h"
#include "multiarray.h"

#include "compute_volume_rlhs.h"
#include "intrusive.h"
#include "multiarray_operator.h"
#include "simulation.h"

// Static function declarations ************************************************************************************* //

/// \brief `complex` version of \ref Flux_Ref.
struct Flux_Ref_c {
	const struct const_Multiarray_c* fr; ///< See brief.
};

/** \brief See \ref compute_volume_rlhs_dg.c.
 *  \return Standard. */
static const struct const_Multiarray_c* constructor_sol_vc
	(struct Solver_Volume* s_vol, ///< See brief.
	 const struct Simulation* sim ///< See brief.
	);

/// \brief See \ref compute_volume_rlhs_dg.c.
static void destructor_sol_vc
	(const struct const_Multiarray_c* sol_vc ///< See brief.
	);

/** \brief `complex` version of \ref constructor_Flux_Ref.
 *  \return See brief. */
static struct Flux_Ref_c* constructor_Flux_Ref_c
	(const struct const_Multiarray_d* m, ///< See brief.
	 const struct Flux_c* flux           ///< See brief.
	);

/// \brief `complex` version of \ref destructor_Flux_Ref.
static void destructor_Flux_Ref_c
	(struct Flux_Ref_c* flux_ref ///< See brief.
	);

/// \brief `complex` version of \ref compute_rhs_v_dg.
static void compute_rhs_v_dg_c
	(const struct Flux_Ref_c* flux_r, ///< See brief.
	 struct Solver_Volume* s_vol,     ///< See brief.
	 const struct Simulation* sim     ///< See brief.
	);

// Interface functions ********************************************************************************************** //

void compute_volume_rhs_dg_c (const struct Simulation* sim, struct Intrusive_List* volumes)
{
	assert(sim->volumes->name == IL_VOLUME_SOLVER_DG_COMPLEX);
	assert(sim->elements->name == IL_ELEMENT_SOLVER_DG);

	struct Flux_Input_c* flux_i = constructor_Flux_Input_c(sim); // destructed

	for (struct Intrusive_Link* curr = volumes->first; curr; curr = curr->next) {
		struct Solver_Volume* s_vol = (struct Solver_Volume*) curr;

		flux_i->s   = constructor_sol_vc(s_vol,sim);
		flux_i->g   = NULL;
		flux_i->xyz = NULL;

		struct Flux_c* flux = constructor_Flux_c(flux_i);
		destructor_sol_vc(flux_i->s);

		struct Flux_Ref_c* flux_r = constructor_Flux_Ref_c(s_vol->metrics_vc,flux);

		destructor_Flux_c(flux);
		compute_rhs_v_dg_c(flux_r,s_vol,sim);
		destructor_Flux_Ref_c(flux_r);
	}
	destructor_Flux_Input_c(flux_i);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

/** \brief `complex` version of \ref constructor_flux_ref.
 *  \return See brief. */
static const struct const_Multiarray_c* constructor_flux_ref_c
	(const struct const_Multiarray_d* m, ///< See brief.
	 const struct const_Multiarray_c* f  ///< See brief.
	);

static const struct const_Multiarray_c* constructor_sol_vc (struct Solver_Volume* s_vol, const struct Simulation* sim)
{
	UNUSED(sim);
	const struct Operator* cv0_vs_vc = get_operator__cv0_vs_vc(s_vol);

	struct Complex_DG_Solver_Volume* dg_s_vol = (struct Complex_DG_Solver_Volume*) s_vol;
	const struct const_Multiarray_c* s_coef = (const struct const_Multiarray_c*) dg_s_vol->sol_coef;

	return constructor_mm_NN1_Operator_const_Multiarray_c(cv0_vs_vc,s_coef,'C','d',s_coef->order,NULL);
}

static void destructor_sol_vc (const struct const_Multiarray_c* sol_vc)
{
	destructor_const_Multiarray_c(sol_vc);
}

static struct Flux_Ref_c* constructor_Flux_Ref_c (const struct const_Multiarray_d* m, const struct Flux_c* flux)
{
	assert(flux->f != NULL);
	assert(m->extents[0] == flux->f->extents[0]);

	struct Flux_Ref_c* flux_r = calloc(1,sizeof *flux_r); // returned

	flux_r->fr = constructor_flux_ref_c(m,flux->f);

	return flux_r;
}

static void destructor_Flux_Ref_c (struct Flux_Ref_c* flux_ref)
{
	destructor_const_Multiarray_c(flux_ref->fr);
	free(flux_ref);
}

static void compute_rhs_v_dg_c
	(const struct Flux_Ref_c* flux_r, struct Solver_Volume* s_vol, const struct Simulation* sim)
{
	const struct Multiarray_Operator tw1_vt_vc = get_operator__tw1_vt_vc(s_vol);

	struct Complex_DG_Solver_Volume* dg_s_volume = (struct Complex_DG_Solver_Volume*) s_vol;
	const ptrdiff_t d = sim->d;
	for (ptrdiff_t dim = 0; dim < d; ++dim)
		mm_NNC_Operator_Multiarray_c(1.0,1.0,tw1_vt_vc.data[dim],flux_r->fr,dg_s_volume->rhs,'d',2,&dim,NULL);
}

// Level 1 ********************************************************************************************************** //

static const struct const_Multiarray_c* constructor_flux_ref_c
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
}
