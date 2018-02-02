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
 *  \brief Provides the templated Navier-Stokes flux functions.
 */

#include <assert.h>
#include <stddef.h>

#include "macros.h"
#include "definitions_core.h"
#include "definitions_physics.h"

#include "def_templates_multiarray.h"

#include "def_templates_flux.h"

// Static function declarations ************************************************************************************* //

#define NEQ  NEQ_EULER  ///< Number of equations.
#define NVAR NVAR_EULER ///< Number of variables.

/// \brief Container for common data used to compute the fluxes and their Jacobians/Hessians.
struct Flux_Data_Navier_Stokes {
	Type rho, ///< The density.
	     E,   ///< The total energy.
	     mu,  ///< The viscosity.
	     Pr;  ///< The Prandtl Number.

	const Type* uvw; ///< The velocity variables.

	const Type*const* tau; ///< Viscous stress tensor.
};

/** \brief Pointer to functions computing the required Navier-Stokes flux members.
 *
 *  \param flux_data  \ref Flux_Data_Navier_Stokes.
 *  \param f_ptr       Pointers to the flux members.
 *  \param dfds_ptr    Pointers to the flux Jacobian members (wrt the solution).
 *  \param dfdg_ptr    Pointers to the flux Jacobian members (wrt the solution gradient).
 *  \todo  d2fds2_ptr  Pointers to the flux Hessian members (wrt the solution).
 *  \todo  d2fdsdg_ptr Pointers to the flux Hessian members (wrt the solution and the gradient).
 *  \todo  d2fdg2_ptr  Pointers to the flux Hessian members (wrt the gradient).
 */
typedef void (*compute_Flux_Navier_Stokes_fptr)
	(const struct Flux_Data_Navier_Stokes*const flux_data,
	 Type*const f_ptr[DIM*NEQ],
	 Type*const dfds_ptr[DIM*NEQ*NVAR],
	 Type*const dfdg_ptr[DIM*NEQ*NVAR*DIM]
	);

/** \brief Return the pointer to the appropriate \ref compute_Flux_Navier_Stokes_fptr specialization based on the
 *         required members.
 *  \return See brief. */
static compute_Flux_Navier_Stokes_fptr get_compute_Flux_Navier_Stokes_fptr
	(const bool*const c_m ///< \ref Flux_Input_T::compute_member.
	);

// Interface functions ********************************************************************************************** //

void compute_Flux_T_navier_stokes (const struct Flux_Input_T* flux_i, struct mutable_Flux_T* flux)
{
	const struct const_Multiarray_T*const s = flux_i->s;
	const Type*const rho_p      = get_col_const_Multiarray_T(0,s),
	          *const rhouvw_p[] = ARRAY_DIM(get_col_const_Multiarray_T(1,s),
	                                        get_col_const_Multiarray_T(2,s),
	                                        get_col_const_Multiarray_T(3,s)),
	          *const E_p        = get_col_const_Multiarray_T(NVAR-1,s);

	const bool* c_m = flux_i->compute_member;

/// \todo Set `compute_flux_navier_stokes_n ` as a member of \ref Flux_Input_T and initialize in the constructor.
	compute_Flux_Navier_Stokes_fptr compute_flux_navier_stokes_n = get_compute_Flux_Navier_Stokes_fptr(c_m);

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

	Type* dfdg_ptr[DIM*NEQ*NVAR*DIM] = { NULL };
	if (c_m[2]) {
		struct Multiarray_T*const dfdg = flux->df_dg;
		for (int dg = 0; dg < DIM; ++dg) {
		for (int vr = 0; vr < NVAR; ++vr)  {
		for (int eq = 0; eq < NEQ; ++eq)  {
		for (int d = 0; d < DIM; ++d) {
			const int ind = d+DIM*(eq+NEQ*(vr+NVAR*(dg)));
			dfdg_ptr[ind] = get_col_Multiarray_T(ind,dfdg);
		}}}}
	}

	assert(c_m[3] == false); // ADD_SUPPORT;
	assert(c_m[4] == false); // ADD_SUPPORT;
	assert(c_m[5] == false); // ADD_SUPPORT;

	const ptrdiff_t n_n = s->extents[0];
	for (ptrdiff_t n = 0; n < n_n; ++n) {
		const Type rho      = rho_p[n],
		           rho_inv  = 1.0/rho,
		           rhouvw[] = ARRAY_DIM(rhouvw_p[0][n],rhouvw_p[1][n],rhouvw_p[2][n]),
		           uvw[]    = ARRAY_DIM(rho_inv*rhouvw[0],rho_inv*rhouvw[1],rho_inv*rhouvw[2]),
		           E        = E_p[n];
EXIT_ADD_SUPPORT; // Compute tau.

		struct Flux_Data_Navier_Stokes flux_data =
			{ .rho = rho,
			  .uvw = uvw,
			  .E   = E,
			};
		compute_flux_navier_stokes_n(&flux_data,f_ptr,dfds_ptr,dfdg_ptr);
	}
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

/// \brief Version of \ref compute_Flux_Navier_Stokes_fptr computing only the flux.
static void compute_Flux_Navier_Stokes_100
	(const struct Flux_Data_Navier_Stokes*const flux_data, ///< See brief.
	 Type*const f_ptr[DIM*NEQ],                            ///< See brief.
	 Type*const dfds_ptr[DIM*NEQ*NVAR],                    ///< See brief.
	 Type*const dfdg_ptr[DIM*NEQ*NVAR*DIM]                 ///< See brief.
	);

/// \brief Version of \ref compute_Flux_Navier_Stokes_fptr computing the flux and Jacobians.
static void compute_Flux_Navier_Stokes_111
	(const struct Flux_Data_Navier_Stokes*const flux_data, ///< See brief.
	 Type*const f_ptr[DIM*NEQ],                            ///< See brief.
	 Type*const dfds_ptr[DIM*NEQ*NVAR],                    ///< See brief.
	 Type*const dfdg_ptr[DIM*NEQ*NVAR*DIM]                 ///< See brief.
	);

static compute_Flux_Navier_Stokes_fptr get_compute_Flux_Navier_Stokes_fptr (const bool*const c_m)
{
	assert(c_m[0]);
	if (c_m[5] || c_m[4] || c_m[3]) {
		EXIT_ADD_SUPPORT;
	} else if (c_m[2]) {
		assert(c_m[1]);
		return compute_Flux_Navier_Stokes_111;
	} else {
		return compute_Flux_Navier_Stokes_100;
	}
}

// Level 1 ********************************************************************************************************** //

/// \brief Compute the Navier-Stokes fluxes for the input nodal values.
static void compute_Flux_navier_stokes_0
	(const struct Flux_Data_Navier_Stokes*const flux_data, ///< \ref Flux_Data_Navier_Stokes.
	 Type*const f_ptr[DIM*NEQ]                             ///< Pointers to the flux data.
	);

/// \brief Compute the Navier-Stokes flux Jacobians wrt the solution for the input nodal values.
static void compute_Flux_navier_stokes_1s
	(const struct Flux_Data_Navier_Stokes*const flux_data, ///< \ref Flux_Data_Navier_Stokes.
	 Type*const dfds_ptr[DIM*NEQ*NVAR]                     ///< Pointers to the flux Jacobian data.
	);

/// \brief Compute the Navier-Stokes flux Jacobians wrt the solution gradients for the input nodal values.
static void compute_Flux_navier_stokes_1g
	(const struct Flux_Data_Navier_Stokes*const flux_data, ///< \ref Flux_Data_Navier_Stokes.
	 Type*const dfdg_ptr[DIM*NEQ*NVAR*DIM]                 ///< Pointers to the flux Jacobian data.
	);

static void compute_Flux_Navier_Stokes_100
	(const struct Flux_Data_Navier_Stokes*const flux_data, Type*const f_ptr[DIM*NEQ],
	 Type*const dfds_ptr[DIM*NEQ*NVAR], Type*const dfdg_ptr[DIM*NEQ*NVAR*DIM])
{
	compute_Flux_navier_stokes_0(flux_data,f_ptr);
	UNUSED(dfds_ptr);
	UNUSED(dfdg_ptr);

	increment_pointers_T(DIM*NEQ,(const Type**)f_ptr);
}

static void compute_Flux_Navier_Stokes_111
	(const struct Flux_Data_Navier_Stokes*const flux_data, Type*const f_ptr[DIM*NEQ],
	 Type*const dfds_ptr[DIM*NEQ*NVAR], Type*const dfdg_ptr[DIM*NEQ*NVAR*DIM])
{
	compute_Flux_navier_stokes_0(flux_data,f_ptr);
	compute_Flux_navier_stokes_1s(flux_data,dfds_ptr);
	compute_Flux_navier_stokes_1g(flux_data,dfdg_ptr);

	increment_pointers_T(DIM*NEQ,         (const Type**)f_ptr);
	increment_pointers_T(DIM*NEQ*NVAR,    (const Type**)dfds_ptr);
	increment_pointers_T(DIM*NEQ*NVAR*DIM,(const Type**)dfdg_ptr);
}

// Level 2 ********************************************************************************************************** //

static void compute_Flux_navier_stokes_0
	(const struct Flux_Data_Navier_Stokes*const flux_data, Type*const f_ptr[DIM*NEQ])
{
	const Type*const*const tau = flux_data->tau;

	IF_DIM_GE_1( const Type u = flux_data->uvw[0]; )
	IF_DIM_GE_2( const Type v = flux_data->uvw[1]; )
	IF_DIM_GE_3( const Type w = flux_data->uvw[2]; )

	const Type mu = flux_data->mu,
	           Pr = flux_data->Pr;

	int ind = 0;

	// f[:,0]
	IF_DIM_GE_1( *f_ptr[ind++] += 0.0 );
	IF_DIM_GE_2( *f_ptr[ind++] += 0.0 );
	IF_DIM_GE_3( *f_ptr[ind++] += 0.0 );

	// f[:,1]
	IF_DIM_GE_1( *f_ptr[ind++] += -tau[0][0] );
	IF_DIM_GE_2( *f_ptr[ind++] += -tau[0][1] );
	IF_DIM_GE_3( *f_ptr[ind++] += -tau[0][2] );

	// f[:,2]
	IF_DIM_GE_2( *f_ptr[ind++] += -tau[1][0] );
	IF_DIM_GE_2( *f_ptr[ind++] += -tau[1][1] );
	IF_DIM_GE_3( *f_ptr[ind++] += -tau[1][2] );

	// f[:,3]
	IF_DIM_GE_3( *f_ptr[ind++] += -tau[2][0] );
	IF_DIM_GE_3( *f_ptr[ind++] += -tau[2][1] );
	IF_DIM_GE_3( *f_ptr[ind++] += -tau[2][2] );

const Type dTs[] = ARRAY_DIM( 0.0, 0.0, 0.0 );
EXIT_ADD_SUPPORT;
	// f[:,4]
	IF_DIM_GE_1( *f_ptr[ind++] += - mu*GAMMA/Pr*dTs[0] - u*tau[0][0] )
	IF_DIM_GE_2(                                       - v*tau[0][1] )
	IF_DIM_GE_3(                                       - w*tau[0][2] );
	IF_DIM_GE_2( *f_ptr[ind++] += - mu*GAMMA/Pr*dTs[1] - u*tau[1][0] )
	IF_DIM_GE_2(                                       - v*tau[1][1] )
	IF_DIM_GE_3(                                       - w*tau[1][2] );
	IF_DIM_GE_3( *f_ptr[ind++] += - mu*GAMMA/Pr*dTs[2] - u*tau[2][0] )
	IF_DIM_GE_3(                                       - v*tau[2][1] )
	IF_DIM_GE_3(                                       - w*tau[2][2] );
}

static void compute_Flux_navier_stokes_1s
	(const struct Flux_Data_Navier_Stokes*const flux_data, Type*const dfds_ptr[DIM*NEQ*NVAR])
{
	EXIT_ADD_SUPPORT; UNUSED(flux_data); UNUSED(dfds_ptr);
}

static void compute_Flux_navier_stokes_1g
	(const struct Flux_Data_Navier_Stokes*const flux_data, Type*const dfdg_ptr[DIM*NEQ*NVAR*DIM])
{
	EXIT_ADD_SUPPORT; UNUSED(flux_data); UNUSED(dfdg_ptr);
}
