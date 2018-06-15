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

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>

#include "macros.h"


#include "def_templates_compute_rlhs.h"
#include "def_templates_multiarray.h"
#include "def_templates_flux.h"
#include "def_templates_math_functions.h"

// Static function declarations ************************************************************************************* //

/** \brief Constructor for a \ref const_Multiarray_T\* of reference flux from physical flux.
 *  \return See brief. */
static const struct const_Multiarray_T* constructor_flux_ref_piece_T
	(const struct const_Multiarray_T*const m, ///< The metric terms.
	 const struct const_Multiarray_T*const f  ///< The physical flux data.
		);

// Interface functions ********************************************************************************************** //

struct Flux_Ref_T* constructor_Flux_Ref_T (const struct const_Multiarray_T*const m, const struct Flux_T*const flux)
{
	assert(flux->f != NULL);
	assert(m->extents[0] == flux->f->extents[0]);

	struct Flux_Ref_T*const flux_r = calloc(1,sizeof *flux_r); // returned

	flux_r->fr       = ( flux->f       ? constructor_flux_ref_piece_T(m,flux->f)       : NULL );
	flux_r->dfr_ds   = ( flux->df_ds   ? constructor_flux_ref_piece_T(m,flux->df_ds)   : NULL );
	flux_r->dfr_dg   = ( flux->df_dg   ? constructor_flux_ref_piece_T(m,flux->df_dg)   : NULL );
	flux_r->d2fr_ds2 = ( flux->d2f_ds2 ? constructor_flux_ref_piece_T(m,flux->d2f_ds2) : NULL );

	return flux_r;
}

void destructor_Flux_Ref_T (struct Flux_Ref_T*const flux_ref)
{
	destructor_conditional_const_Multiarray_T(flux_ref->fr);
	destructor_conditional_const_Multiarray_T(flux_ref->dfr_ds);
	destructor_conditional_const_Multiarray_T(flux_ref->dfr_dg);
	destructor_conditional_const_Multiarray_T(flux_ref->d2fr_ds2);
	free(flux_ref);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

static const struct const_Multiarray_T* constructor_flux_ref_piece_T
	(const struct const_Multiarray_T*const m, const struct const_Multiarray_T*const f)
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

	struct Multiarray_T* fr = constructor_zero_Multiarray_T('C',order,extents); // returned

	const int n_n   = (int)extents[0];
	const int n_col = (int)compute_size(order,extents)/(n_n*DIM);
	assert(extents[order-1] == DIM);

	int ind_f = 0;
	for (int col = 0; col < n_col; ++col) {
		for (int dim0 = 0; dim0 < DIM; ++dim0) {
			const int ind_fr = (ind_f+dim0*n_col);
			for (int dim1 = 0; dim1 < DIM; ++dim1) {
				const int ind_m  = dim0*DIM+dim1,
				          ind_fp = (ind_f*DIM)+dim1;
				z_yxpz_T(n_n,
				         get_col_const_Multiarray_T(ind_m,m),
				         get_col_const_Multiarray_T(ind_fp,f),
				         get_col_Multiarray_T(ind_fr,fr));
			}
		}
		++ind_f;
	}

	return (const struct const_Multiarray_T*) fr;
}

#include "undef_templates_compute_rlhs.h"
#include "undef_templates_multiarray.h"
#include "undef_templates_flux.h"
#include "undef_templates_math_functions.h"
