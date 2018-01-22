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
 *  \brief Provides the templated diffusion flux functions.
 */

#include <assert.h>
#include <stddef.h>

#include "macros.h"
#include "definitions_core.h"

#include "def_templates_multiarray.h"

#include "def_templates_flux.h"

// Static function declarations ************************************************************************************* //

#define NEQ  1 ///< Number of equations.
#define NVAR 1 ///< Number of variables.

/// \brief Container for common data used to compute the fluxes and their Jacobians.
struct Flux_Data_Diffusion {
	const Type* g; ///< Gradient variable.
};

/** \brief Pointer to functions computing the required Diffusion flux members.
 *
 *  \param flux_data  \ref Flux_Data_Diffusion.
 *  \param f_ptr      Pointers to the flux members.
 *  \param dfdg_ptr   Pointers to the flux Jacobian members.
 */
typedef void (*compute_Flux_Diffusion_fptr)
	(const struct Flux_Data_Diffusion*const flux_data,
	 Type*const f_ptr[DIM*NEQ],
	 Type*const dfdg_ptr[DIM*DIM*NEQ*NVAR]
	);

/** \brief Return the pointer to the appropriate \ref compute_Flux_Diffusion_fptr specialization based on the required
 *         members.
 *  \return See brief. */
static compute_Flux_Diffusion_fptr get_compute_Flux_Diffusion_fptr
	(const bool*const c_m ///< \ref Flux_Input_T::compute_member.
	);

// Interface functions ********************************************************************************************** //

void compute_Flux_T_diffusion (const struct Flux_Input_T* flux_i, struct mutable_Flux_T* flux)
{
	const struct const_Multiarray_T*const s = flux_i->s;
	const struct const_Multiarray_T*const g = flux_i->g;

	const ptrdiff_t size_s = compute_size(s->order,s->extents);

	const Type*const g_p[] = ARRAY_DIM( get_col_const_Multiarray_T(0*size_s,g),
	                                    get_col_const_Multiarray_T(1*size_s,g),
	                                    get_col_const_Multiarray_T(2*size_s,g) );

	const bool* c_m = flux_i->compute_member;

	compute_Flux_Diffusion_fptr compute_flux_diffusion_n = get_compute_Flux_Diffusion_fptr(c_m);

	assert(c_m[0]);
	Type* f_ptr[DIM*NEQ] = { NULL };
	struct Multiarray_T*const f = flux->f;
	for (int eq = 0; eq < NEQ; ++eq)  {
	for (int d = 0; d < DIM; ++d) {
		const int ind = d+DIM*(eq);
		f_ptr[ind] = get_col_Multiarray_T(ind,f);
	}}

	assert(c_m[1] == false); // 2nd order only.

	Type* dfdg_ptr[DIM*NEQ*NVAR*DIM] = { NULL };
	if (c_m[2]) {
		struct Multiarray_T*const dfdg = flux->df_dg;
		for (int dd = 0; dd < DIM; ++dd) {
		for (int vr = 0; vr < NVAR; ++vr)  {
		for (int eq = 0; eq < NEQ; ++eq)  {
		for (int d = 0; d < DIM; ++d) {
			const int ind = d+DIM*(eq+NEQ*(vr+NVAR*(dd)));
			dfdg_ptr[ind] = get_col_Multiarray_T(ind,dfdg);
		}}}}
	}

	assert(c_m[3] == false); // linear.
	assert(c_m[4] == false); // linear.
	assert(c_m[5] == false); // linear.

	const ptrdiff_t n_n = s->extents[0];
	for (ptrdiff_t n = 0; n < n_n; ++n) {
		Type g_n[] = ARRAY_DIM ( g_p[0][n], g_p[1][n], g_p[2][n] );
		struct Flux_Data_Diffusion flux_data =
			{ .g = g_n,
			};
		compute_flux_diffusion_n(&flux_data,f_ptr,dfdg_ptr);
	}
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

/// \brief Version of \ref compute_Flux_Diffusion_fptr computing only the flux.
static void compute_Flux_Diffusion_100
	(const struct Flux_Data_Diffusion*const flux_data, ///< See brief.
	 Type*const f_ptr[DIM*NEQ],                        ///< See brief.
	 Type*const dfdg_ptr[DIM*NEQ*NVAR*DIM]             ///< See brief.
	);

/// \brief Version of \ref compute_Flux_Diffusion_fptr computing the flux and Jacobian.
static void compute_Flux_Diffusion_101
	(const struct Flux_Data_Diffusion*const flux_data, ///< See brief.
	 Type*const f_ptr[DIM*NEQ],                        ///< See brief.
	 Type*const dfdg_ptr[DIM*NEQ*NVAR*DIM]             ///< See brief.
	);

static compute_Flux_Diffusion_fptr get_compute_Flux_Diffusion_fptr (const bool*const c_m)
{
	assert(c_m[0]);
	if (c_m[2])
		return compute_Flux_Diffusion_101;
	else
		return compute_Flux_Diffusion_100;
}

// Level 1 ********************************************************************************************************** //

/** \brief Compute the Diffusion fluxes for the input nodal values.
 *  \warning  Negated flux contributions are returned such that the identical treatment to inviscid flux contributions
 *            is applicable. */
static void compute_Flux_diffusion_0
	(const struct Flux_Data_Diffusion*const flux_data, ///< \ref Flux_Data_Diffusion.
	 Type*const f_ptr[DIM*NEQ]                         ///< Pointers to the flux data.
	);

/** \brief Compute the Diffusion flux Jacobians for the input nodal values.
 *  \warning  Negated flux contributions are returned such that the identical treatment to inviscid flux contributions
 *            is applicable. */
static void compute_Flux_diffusion_2
	(const struct Flux_Data_Diffusion*const flux_data, ///< \ref Flux_Data_Diffusion.
	 Type*const dfdg_ptr[DIM*NEQ*NVAR*DIM]             ///< Pointers to the flux Jacobian data.
	);

static void compute_Flux_Diffusion_100
	(const struct Flux_Data_Diffusion*const flux_data, Type*const f_ptr[DIM*NEQ],
	 Type*const dfdg_ptr[DIM*NEQ*NVAR*DIM])
{
	compute_Flux_diffusion_0(flux_data,f_ptr);
	UNUSED(dfdg_ptr);

	increment_pointers_T(DIM*NEQ,(const Type**)f_ptr);
}

static void compute_Flux_Diffusion_101
	(const struct Flux_Data_Diffusion*const flux_data, Type*const f_ptr[DIM*NEQ],
	 Type*const dfdg_ptr[DIM*NEQ*NVAR*DIM])
{
	compute_Flux_diffusion_0(flux_data,f_ptr);
	compute_Flux_diffusion_2(flux_data,dfdg_ptr);

	increment_pointers_T(DIM*NEQ,         (const Type**)f_ptr);
	increment_pointers_T(DIM*NEQ*NVAR*DIM,(const Type**)dfdg_ptr);
}

// Level 2 ********************************************************************************************************** //

static void compute_Flux_diffusion_0 (const struct Flux_Data_Diffusion*const flux_data, Type*const f_ptr[DIM*NEQ])
{
	const Type* g = flux_data->g;

	int ind = 0;

	// f[:,0]
	IF_DIM_GE_1( *f_ptr[ind++] += -g[0] );
	IF_DIM_GE_2( *f_ptr[ind++] += -g[1] );
	IF_DIM_GE_3( *f_ptr[ind++] += -g[2] );
}

static void compute_Flux_diffusion_2
	(const struct Flux_Data_Diffusion*const flux_data, Type*const dfdg_ptr[DIM*NEQ*NVAR*DIM])
{
	UNUSED(flux_data);

	int ind = 0;

	// dfdg[:,:,0]
	// dfdg[:,0,0]
	IF_DIM_GE_1( *dfdg_ptr[ind++] += -1.0 );
	IF_DIM_GE_2( *dfdg_ptr[ind++] +=  0.0 );
	IF_DIM_GE_3( *dfdg_ptr[ind++] +=  0.0 );

	// dfdg[:,0,1]
	IF_DIM_GE_1( *dfdg_ptr[ind++] +=  0.0 );
	IF_DIM_GE_2( *dfdg_ptr[ind++] += -1.0 );
	IF_DIM_GE_3( *dfdg_ptr[ind++] +=  0.0 );

	// dfdg[:,0,2]
	IF_DIM_GE_1( *dfdg_ptr[ind++] +=  0.0 );
	IF_DIM_GE_2( *dfdg_ptr[ind++] +=  0.0 );
	IF_DIM_GE_3( *dfdg_ptr[ind++] += -1.0 );
}
