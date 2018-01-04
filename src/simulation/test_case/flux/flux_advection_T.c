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
 *  \brief Provides the templated linear advection flux functions.
 */

#include <assert.h>
#include <stddef.h>

#include "macros.h"

#include "def_templates_multiarray.h"

#include "def_templates_flux.h"

// Static function declarations ************************************************************************************* //

#define NEQ  1 ///< Number of equations.
#define NVAR 1 ///< Number of variables.

/// \brief Container for common data used to compute the fluxes and their Jacobians.
struct Flux_Data_Advection {
	const Type u; ///< State variable.

	const Real* b_adv; ///< Advection velocity vector.
};

/** \brief Pointer to functions computing the required Advection flux members.
 *
 *  \param flux_data  \ref Flux_Data_Advection.
 *  \param f_ptr      Pointers to the flux members.
 *  \param dfds_ptr   Pointers to the flux Jacobian members.
 */
typedef void (*compute_Flux_Advection_fptr)
	(const struct Flux_Data_Advection*const flux_data,
	 Type*const f_ptr[DIM*NEQ],
	 Type*const dfds_ptr[DIM*NEQ*NVAR]
	);

/** \brief Return the pointer to the appropriate \ref compute_Flux_Advection_fptr specialization based on the required
 *         members.
 *  \return See brief. */
static compute_Flux_Advection_fptr get_compute_Flux_Advection_fptr
	(const bool*const c_m ///< \ref Flux_Input_T::compute_member.
	);

// Interface functions ********************************************************************************************** //

void compute_Flux_T_advection (const struct Flux_Input_T* flux_i, struct mutable_Flux_T* flux)
{
	static bool need_input = true;
	static struct Sol_Data__Advection sol_data;
	if (need_input) {
		need_input = false;
		read_data_advection(flux_i->input_path,&sol_data);
	}

	const struct const_Multiarray_T*const s = flux_i->s;
	const Type*const u_p = get_col_const_Multiarray_T(0,s);

	const bool* c_m = flux_i->compute_member;

	compute_Flux_Advection_fptr compute_flux_advection_n = get_compute_Flux_Advection_fptr(c_m);

	assert(c_m[0]);
	Type* f_ptr[DIM*NEQ] = { NULL };
	struct Multiarray_T*const f = flux->f;
	for (int eq = 0; eq < NEQ; ++eq)  {
	for (int d = 0; d < DIM; ++d) {
		const int ind = d+DIM*(eq);
		f_ptr[ind] = get_col_Multiarray_T(ind,f);
	}}

	Type* dfds_ptr[DIM*NEQ*NVAR] = { NULL };
	if (c_m[1]) {
		struct Multiarray_T*const dfds = flux->df_ds;
		for (int vr = 0; vr < NVAR; ++vr)  {
		for (int eq = 0; eq < NEQ; ++eq)  {
		for (int d = 0; d < DIM; ++d) {
			const int ind = d+DIM*(eq+NEQ*(vr));
			dfds_ptr[ind] = get_col_Multiarray_T(ind,dfds);
		}}}
	}

	assert(c_m[2] == false); // 1st order.
	assert(c_m[3] == false); // linear.
	assert(c_m[4] == false); // 1st order.
	assert(c_m[5] == false); // 1st order.

	const ptrdiff_t n_n = s->extents[0];
	for (ptrdiff_t n = 0; n < n_n; ++n) {
		struct Flux_Data_Advection flux_data =
			{ .u     = u_p[n],
			  .b_adv = sol_data.b_adv,
			};
		compute_flux_advection_n(&flux_data,f_ptr,dfds_ptr);
	}
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

/// \brief Version of \ref compute_Flux_Advection_fptr computing only the flux.
static void compute_Flux_Advection_100
	(const struct Flux_Data_Advection*const flux_data, ///< See brief.
	 Type*const f_ptr[DIM*NEQ],                        ///< See brief.
	 Type*const dfds_ptr[DIM*NEQ*NVAR]                 ///< See brief.
	);

/// \brief Version of \ref compute_Flux_Advection_fptr computing the flux and Jacobian.
static void compute_Flux_Advection_110
	(const struct Flux_Data_Advection*const flux_data, ///< See brief.
	 Type*const f_ptr[DIM*NEQ],                        ///< See brief.
	 Type*const dfds_ptr[DIM*NEQ*NVAR]                 ///< See brief.
	);

static compute_Flux_Advection_fptr get_compute_Flux_Advection_fptr (const bool*const c_m)
{
	assert(c_m[0]);
	if (c_m[1])
		return compute_Flux_Advection_110;
	else
		return compute_Flux_Advection_100;
}

// Level 1 ********************************************************************************************************** //

/// \brief Compute the Advection fluxes for the input nodal values.
static void compute_Flux_advection_0
	(const struct Flux_Data_Advection*const flux_data, ///< \ref Flux_Data_Advection.
	 Type*const f_ptr[DIM*NEQ]                         ///< Pointers to the flux data.
	);

/// \brief Compute the Advection flux Jacobians for the input nodal values.
static void compute_Flux_advection_1
	(const struct Flux_Data_Advection*const flux_data, ///< \ref Flux_Data_Advection.
	 Type*const dfds_ptr[DIM*NEQ*NVAR]                 ///< Pointers to the flux Jacobian data.
	);

static void compute_Flux_Advection_100
	(const struct Flux_Data_Advection*const flux_data, Type*const f_ptr[DIM*NEQ], Type*const dfds_ptr[DIM*NEQ*NVAR])
{
	compute_Flux_advection_0(flux_data,f_ptr);
	UNUSED(dfds_ptr);

	increment_pointers_T(DIM*NEQ,(const Type**)f_ptr);
}

static void compute_Flux_Advection_110
	(const struct Flux_Data_Advection*const flux_data, Type*const f_ptr[DIM*NEQ], Type*const dfds_ptr[DIM*NEQ*NVAR])
{
	compute_Flux_advection_0(flux_data,f_ptr);
	compute_Flux_advection_1(flux_data,dfds_ptr);

	increment_pointers_T(DIM*NEQ,     (const Type**)f_ptr);
	increment_pointers_T(DIM*NEQ*NVAR,(const Type**)dfds_ptr);
}

// Level 2 ********************************************************************************************************** //

static void compute_Flux_advection_0 (const struct Flux_Data_Advection*const flux_data, Type*const f_ptr[DIM*NEQ])
{
	const Type u = flux_data->u;

	const Real*const b_adv = flux_data->b_adv;

	int ind = 0;

	// f[:,0]
	IF_DIM_GE_1( *f_ptr[ind++] += b_adv[0]*u );
	IF_DIM_GE_2( *f_ptr[ind++] += b_adv[1]*u );
	IF_DIM_GE_3( *f_ptr[ind++] += b_adv[2]*u );
}

static void compute_Flux_advection_1
	(const struct Flux_Data_Advection*const flux_data, Type*const dfds_ptr[DIM*NEQ*NVAR])
{
	const Real*const b_adv = flux_data->b_adv;

	int ind = 0;

	// dfds[:,:,0]
	// dfds[:,0,0]
	IF_DIM_GE_1( *dfds_ptr[ind++] += b_adv[0] );
	IF_DIM_GE_2( *dfds_ptr[ind++] += b_adv[1] );
	IF_DIM_GE_3( *dfds_ptr[ind++] += b_adv[2] );
}
