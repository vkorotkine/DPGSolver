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

#include "test_complex_compute_volume_rhs.h"

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>

#include "macros.h"
#include "definitions_intrusive.h"

#include "test_complex_compute_volume_rhs_dg.h"

#include "matrix.h"
#include "multiarray.h"
#include "vector.h"

#include "volume_solver.h"
#include "element_solver.h"

#include "flux.h"
#include "intrusive.h"
#include "math_functions.h"
#include "multiarray_operator.h"
#include "operator.h"
#include "simulation.h"
#include "test_case.h"

// Static function declarations ************************************************************************************* //

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

// Interface functions ********************************************************************************************** //

struct Flux_Ref_c* constructor_Flux_Ref_vol_c
	(const struct S_Params_Volume_Structor_c* spvs, struct Flux_Input_c* flux_i, const struct Solver_Volume* s_vol,
	 const struct Simulation* sim)
{
	flux_i->s   = spvs->constructor_sol_vc(s_vol,sim);
	flux_i->g   = NULL;
	flux_i->xyz = NULL;

	struct Flux_c* flux = constructor_Flux_c(flux_i);
	destructor_sol_vc(flux_i->s);

	struct Flux_Ref_c* flux_r = constructor_Flux_Ref_c(s_vol->metrics_vc,flux);
	destructor_Flux_c(flux);

	return flux_r;
}

void destructor_Flux_Ref_c (struct Flux_Ref_c* flux_ref)
{
	destructor_const_Multiarray_c(flux_ref->fr);
	free(flux_ref);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

/** \brief `complex` version of \ref constructor_flux_ref.
 *  \return See brief. */
static const struct const_Multiarray_c* constructor_flux_ref_c
	(const struct const_Multiarray_d* m, ///< See brief.
	 const struct const_Multiarray_c* f  ///< See brief.
	);

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
