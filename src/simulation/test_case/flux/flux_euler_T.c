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
 *  \brief Provides the templated Euler flux functions.
 */

#include <assert.h>
#include <stddef.h>

#include "macros.h"
#include "definitions_core.h"
#include "definitions_test_case.h"

#include "def_templates_flux.h"

// Static function declarations ************************************************************************************* //

#define NEQ  NEQ_EULER  ///< Number of equations.
#define NVAR NVAR_EULER ///< Number of variables.

/// \brief Container for common data used to compute the fluxes and their Jacobians/Hessians.
struct Flux_Data_Euler {
	const Type* rhouvw, ///< Array of momentum variables.
	          * uvw;    ///< Array of velocity variables.
	const Type rho_inv, ///< Inverse of density.
	           E,       ///< Total energy.
		   p,       ///< Pressure.
		   V2;      ///< Square of the velocity magnitude.
};

/** \brief Pointer to functions computing the required Euler flux members.
 *
 *  \param flux_data  \ref Flux_Data_Euler.
 *  \param f_ptr      Pointers to the flux members.
 *  \param dfds_ptr   Pointers to the flux Jacobian members.
 *  \param d2fds2_ptr Pointers to the flux Hessian members.
 */
typedef void (*compute_Flux_Euler_fptr)
	(const struct Flux_Data_Euler*const flux_data, ///< \ref Flux_Data_Euler.
	 Type*const f_ptr[DIM*NEQ],                    ///< Pointers to the flux members.
	 Type*const dfds_ptr[DIM*NEQ*NVAR],            ///< Pointers to the flux Jacobian members.
	 Type*const d2fds2_ptr[DIM*NEQ*NVAR*NVAR]      ///< Pointers to the flux Hessian members.
	);

/** \brief Return the pointer to the appropriate \ref compute_Flux_Euler_fptr specialization based on the required
 *         members.
 *  \return See brief. */
static compute_Flux_Euler_fptr get_compute_Flux_Euler_fptr
	(const bool*const c_m ///< \ref Flux_Input_T::compute_member.
	);

/** \brief Compute the square of the Velocity magnitude.
 *  \return See brief. */
static Type compute_V2
	(const Type* uvw ///< The array of velocity components.
	);

// Interface functions ********************************************************************************************** //

void compute_Flux_T_euler (const struct Flux_Input_T* flux_i, struct mutable_Flux_T* flux)
{
	const struct const_Multiarray_T*const s = flux_i->s;
	const Type*const rho_p      = get_col_const_Multiarray_T(0,s),
	          *const rhouvw_p[] = ARRAY_DIM(get_col_const_Multiarray_T(1,s),
	                                        get_col_const_Multiarray_T(2,s),
	                                        get_col_const_Multiarray_T(3,s)),
	          *const E_p        = get_col_const_Multiarray_T(NVAR-1,s);

	const bool* c_m = flux_i->compute_member;

	compute_Flux_Euler_fptr compute_flux_euler_n = get_compute_Flux_Euler_fptr(c_m);

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

	Type* d2fds2_ptr[DIM*NEQ*NVAR*NVAR] = { NULL };
	if (c_m[3]) {
		struct Multiarray_T*const d2fds2 = flux->d2f_ds2;
		for (int vr2 = 0; vr2 < NVAR; ++vr2)  {
		for (int vr = 0; vr < NVAR; ++vr)  {
		for (int eq = 0; eq < NEQ; ++eq)  {
		for (int d = 0; d < DIM; ++d) {
			const int ind = d+DIM*(eq+NEQ*(vr+NVAR*(vr2)));
			d2fds2_ptr[ind] = get_col_Multiarray_T(ind,d2fds2);
		}}}}
	}

	const ptrdiff_t n_n = s->extents[0];
	for (ptrdiff_t n = 0; n < n_n; ++n) {
		const Type rho      = rho_p[n],
		           rho_inv  = 1.0/rho,

		           rhouvw[] = ARRAY_DIM(rhouvw_p[0][n],rhouvw_p[1][n],rhouvw_p[2][n]),
		           uvw[]    = ARRAY_DIM(rho_inv*rhouvw[0],rho_inv*rhouvw[1],rho_inv*rhouvw[2]),

		           V2 = compute_V2(uvw),
		           E = E_p[n],
		           p = GM1*(E-0.5*rho*V2);

		struct Flux_Data_Euler flux_data =
			{ .rhouvw  = rhouvw,
			  .uvw     = uvw,
			  .rho_inv = rho_inv,
			  .E       = E,
			  .p       = p,
			  .V2      = V2,
			};
		compute_flux_euler_n(&flux_data,f_ptr,dfds_ptr,d2fds2_ptr);
	}
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

/// \brief Version of \ref compute_Flux_Euler_fptr computing only the flux.
static void compute_Flux_Euler_100
	(const struct Flux_Data_Euler*const flux_data, ///< See brief.
	 Type*const f_ptr[DIM*NEQ],                    ///< See brief.
	 Type*const dfds_ptr[DIM*NEQ*NVAR],            ///< See brief.
	 Type*const d2fds2_ptr[DIM*NEQ*NVAR*NVAR]      ///< See brief.
	);

/// \brief Version of \ref compute_Flux_Euler_fptr computing the flux and Jacobian.
static void compute_Flux_Euler_110
	(const struct Flux_Data_Euler*const flux_data, ///< See brief.
	 Type*const f_ptr[DIM*NEQ],                    ///< See brief.
	 Type*const dfds_ptr[DIM*NEQ*NVAR],            ///< See brief.
	 Type*const d2fds2_ptr[DIM*NEQ*NVAR*NVAR]      ///< See brief.
	);

/// \brief Version of \ref compute_Flux_Euler_fptr computing the flux, Jacobian and Hessian.
static void compute_Flux_Euler_111
	(const struct Flux_Data_Euler*const flux_data, ///< See brief.
	 Type*const f_ptr[DIM*NEQ],                    ///< See brief.
	 Type*const dfds_ptr[DIM*NEQ*NVAR],            ///< See brief.
	 Type*const d2fds2_ptr[DIM*NEQ*NVAR*NVAR]      ///< See brief.
	);

static compute_Flux_Euler_fptr get_compute_Flux_Euler_fptr (const bool*const c_m)
{
	assert(c_m[0]);
	if (c_m[3]) {
		assert(c_m[1]);
		return compute_Flux_Euler_111;
	} else if (c_m[1]) {
		return compute_Flux_Euler_110;
	} else {
		return compute_Flux_Euler_100;
	}
}

static Type compute_V2 (const Type* uvw)
{
	Type V2 = 0.0;
	for (int d = 0; d < DIM; ++d)
		V2 += uvw[d]*uvw[d];
	return V2;
}

// Level 1 ********************************************************************************************************** //

/// \brief Compute the Euler fluxes for the input nodal values.
static void compute_Flux_euler_0
	(const struct Flux_Data_Euler*const flux_data, ///< \ref Flux_Data_Euler.
	 Type*const f_ptr[DIM*NEQ]                     ///< Pointers to the flux data.
	);

/// \brief Compute the Euler flux Jacobians for the input nodal values.
static void compute_Flux_euler_1
	(const struct Flux_Data_Euler*const flux_data, ///< \ref Flux_Data_Euler.
	 Type*const dfds_ptr[DIM*NEQ*NVAR]             ///< Pointers to the flux Jacobian data.
	);

/// \brief Compute the Euler flux Jacobians for the input nodal values.
static void compute_Flux_euler_2
	(const struct Flux_Data_Euler*const flux_data, ///< \ref Flux_Data_Euler.
	 Type*const d2fds2_ptr[DIM*NEQ*NVAR*NVAR]      ///< Pointers to the flux Hessian data.
	);

/// \brief Increment input number of pointers by one.
static void increment_pointers
	(const int n_ptr,       ///< The number of pointers.
	 const Type**const ptrs ///< Pointer to the array of pointers.
	);

static void compute_Flux_Euler_100
	(const struct Flux_Data_Euler*const flux_data, Type*const f_ptr[DIM*NEQ],
	 Type*const dfds_ptr[DIM*NEQ*NVAR], Type*const d2fds2_ptr[DIM*NEQ*NVAR*NVAR])
{
	compute_Flux_euler_0(flux_data,f_ptr);
	UNUSED(dfds_ptr);
	UNUSED(d2fds2_ptr);

	increment_pointers(DIM*NEQ,(const Type**)f_ptr);
}

static void compute_Flux_Euler_110
	(const struct Flux_Data_Euler*const flux_data, Type*const f_ptr[DIM*NEQ],
	 Type*const dfds_ptr[DIM*NEQ*NVAR], Type*const d2fds2_ptr[DIM*NEQ*NVAR*NVAR])
{
	compute_Flux_euler_0(flux_data,f_ptr);
	compute_Flux_euler_1(flux_data,dfds_ptr);
	UNUSED(d2fds2_ptr);

	increment_pointers(DIM*NEQ,     (const Type**)f_ptr);
	increment_pointers(DIM*NEQ*NVAR,(const Type**)dfds_ptr);
}

static void compute_Flux_Euler_111
	(const struct Flux_Data_Euler*const flux_data, Type*const f_ptr[DIM*NEQ],
	 Type*const dfds_ptr[DIM*NEQ*NVAR], Type*const d2fds2_ptr[DIM*NEQ*NVAR*NVAR])
{
	compute_Flux_euler_0(flux_data,f_ptr);
	compute_Flux_euler_1(flux_data,dfds_ptr);
	compute_Flux_euler_2(flux_data,d2fds2_ptr);

	increment_pointers(DIM*NEQ,          (const Type**)f_ptr);
	increment_pointers(DIM*NEQ*NVAR,     (const Type**)dfds_ptr);
	increment_pointers(DIM*NEQ*NVAR*NVAR,(const Type**)d2fds2_ptr);
}

// Level 2 ********************************************************************************************************** //

static void compute_Flux_euler_0 (const struct Flux_Data_Euler*const flux_data, Type*const f_ptr[DIM*NEQ])
{
	const Type p = flux_data->p,
	           E = flux_data->E;

	const Type*const rhouvw = flux_data->rhouvw,
	          *const uvw    = flux_data->uvw;

	int ind = 0;
	for (int d = 0; d < DIM; ++d)
		*f_ptr[ind++] += rhouvw[d];

	for (int ind_m = 0; ind_m < DIM; ++ind_m) {
	for (int d = 0; d < DIM; ++d) {
		if (d == ind_m)
			*f_ptr[ind++] += rhouvw[ind_m]*uvw[d] + p;
		else
			*f_ptr[ind++] += rhouvw[ind_m]*uvw[d];
	}}

	for (int d = 0; d < DIM; ++d)
		*f_ptr[ind++] += (E+p)*uvw[d];
}

static void compute_Flux_euler_1 (const struct Flux_Data_Euler*const flux_data, Type*const dfds_ptr[DIM*NEQ*NVAR])
{
	const Type rho_inv = flux_data->rho_inv,
	           E       = flux_data->E,
	           p       = flux_data->p,
	           V2      = flux_data->V2;

	const Type*const uvw = flux_data->uvw;

	const Type H     = (E+p)*rho_inv,
	           alpha = 0.5*GM1*V2,
	           beta  = alpha-H;

	int ind = 0;

	// dfds[:,:,0]
	for (int d = 0; d < DIM; ++d)
		*dfds_ptr[ind++] += 0.0;

	for (int ind_m = 0; ind_m < DIM; ++ind_m) {
	for (int d = 0; d < DIM; ++d) {
		if (d == ind_m)
			*dfds_ptr[ind++] += -uvw[ind_m]*uvw[d]+alpha;
		else
			*dfds_ptr[ind++] += -uvw[ind_m]*uvw[d];
	}}

	for (int d = 0; d < DIM; ++d)
		*dfds_ptr[ind++] += uvw[d]*beta;

	// dfds[:,:,1:DIM]
	for (int ind_v = 0; ind_v < DIM; ++ind_v) {
		for (int d = 0; d < DIM; ++d) {
			if (d == ind_v)
				*dfds_ptr[ind++] += 1.0;
			else
				*dfds_ptr[ind++] += 0.0;
		}

		for (int ind_m = 0; ind_m < DIM; ++ind_m) {
		for (int d = 0; d < DIM; ++d) {
			if (ind_m == ind_v) {
				if (ind_m == d)
					*dfds_ptr[ind++] += -GM3*uvw[d];
				else
					*dfds_ptr[ind++] += uvw[d];
			} else {
				if (ind_m == d)
					*dfds_ptr[ind++] += -GM1*uvw[ind_v];
				else if (ind_v == d)
					*dfds_ptr[ind++] += uvw[ind_m];
				else
					*dfds_ptr[ind++] += 0.0;
			}
		}}

		for (int d = 0; d < DIM; ++d) {
			if (d == ind_v)
				*dfds_ptr[ind++] += H-GM1*uvw[ind_v]*uvw[d];
			else
				*dfds_ptr[ind++] += -GM1*uvw[ind_v]*uvw[d];
		}
	}

	// dfds[:,:,DIM+1]
	for (int d = 0; d < DIM; ++d)
		*dfds_ptr[ind++] += 0.0;

	for (int ind_m = 0; ind_m < DIM; ++ind_m) {
	for (int d = 0; d < DIM; ++d) {
		if (d == ind_m)
			*dfds_ptr[ind++] += GM1;
		else
			*dfds_ptr[ind++] += 0.0;
	}}

	for (int d = 0; d < DIM; ++d)
		*dfds_ptr[ind++] += GAMMA*uvw[d];
}

static void compute_Flux_euler_2
	(const struct Flux_Data_Euler*const flux_data, Type*const d2fds2_ptr[DIM*NEQ*NVAR*NVAR])
{
UNUSED(flux_data);
UNUSED(d2fds2_ptr);
printf("Not yet implemented.\n");
}

static void increment_pointers (const int n_ptr, const Type**const ptrs)
{
	for (int i = 0; i < n_ptr; ++i)
		++ptrs[i];
}
