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
	           V2,      ///< Square of the velocity magnitude.
		     H,       ///< Enthalpy.
		     alpha,   ///< dp_drho.
		     beta;    ///< alpha-H.
};

/** \brief Pointer to functions computing the required Euler flux members.
 *
 *  \param flux_data  \ref Flux_Data_Euler.
 *  \param f_ptr      Pointers to the flux members.
 *  \param dfds_ptr   Pointers to the flux Jacobian members.
 *  \param d2fds2_ptr Pointers to the flux Hessian members.
 */
typedef void (*compute_Flux_Euler_fptr)
	(const struct Flux_Data_Euler*const flux_data,
	 Type*const f_ptr[DIM*NEQ],
	 Type*const dfds_ptr[DIM*NEQ*NVAR],
	 Type*const d2fds2_ptr[DIM*NEQ*NVAR*NVAR]
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

/// \todo Set `compute_flux_euler_n ` as a member of \ref Flux_Input_T and initialize in the constructor.
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

	assert(c_m[2] == false); // 1st order.

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

	assert(c_m[4] == false); // 1st order.
	assert(c_m[5] == false); // 1st order.

	const ptrdiff_t n_n = s->extents[0];
	for (ptrdiff_t n = 0; n < n_n; ++n) {
		const Type rho      = rho_p[n],
		           rho_inv  = 1.0/rho,

		           rhouvw[] = ARRAY_DIM(rhouvw_p[0][n],rhouvw_p[1][n],rhouvw_p[2][n]),
		           uvw[]    = ARRAY_DIM(rho_inv*rhouvw[0],rho_inv*rhouvw[1],rho_inv*rhouvw[2]),

		           V2    = compute_V2(uvw),
		           E     = E_p[n],
		           p     = GM1*(E-0.5*rho*V2),
		           H     = (E+p)*rho_inv,
		           alpha = 0.5*GM1*V2;

		struct Flux_Data_Euler flux_data =
			{ .rhouvw  = rhouvw,
			  .uvw     = uvw,
			  .rho_inv = rho_inv,
			  .E       = E,
			  .p       = p,
			  .V2      = V2,
			  .H       = H,
			  .alpha   = alpha,
			  .beta    = alpha-H,
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

/// \brief Compute the Euler flux Hessians for the input nodal values.
static void compute_Flux_euler_2
	(const struct Flux_Data_Euler*const flux_data, ///< \ref Flux_Data_Euler.
	 Type*const d2fds2_ptr[DIM*NEQ*NVAR*NVAR]      ///< Pointers to the flux Hessian data.
	);

static void compute_Flux_Euler_100
	(const struct Flux_Data_Euler*const flux_data, Type*const f_ptr[DIM*NEQ],
	 Type*const dfds_ptr[DIM*NEQ*NVAR], Type*const d2fds2_ptr[DIM*NEQ*NVAR*NVAR])
{
	compute_Flux_euler_0(flux_data,f_ptr);
	UNUSED(dfds_ptr);
	UNUSED(d2fds2_ptr);

	increment_pointers_T(DIM*NEQ,(const Type**)f_ptr);
}

static void compute_Flux_Euler_110
	(const struct Flux_Data_Euler*const flux_data, Type*const f_ptr[DIM*NEQ],
	 Type*const dfds_ptr[DIM*NEQ*NVAR], Type*const d2fds2_ptr[DIM*NEQ*NVAR*NVAR])
{
	compute_Flux_euler_0(flux_data,f_ptr);
	compute_Flux_euler_1(flux_data,dfds_ptr);
	UNUSED(d2fds2_ptr);

	increment_pointers_T(DIM*NEQ,     (const Type**)f_ptr);
	increment_pointers_T(DIM*NEQ*NVAR,(const Type**)dfds_ptr);
}

static void compute_Flux_Euler_111
	(const struct Flux_Data_Euler*const flux_data, Type*const f_ptr[DIM*NEQ],
	 Type*const dfds_ptr[DIM*NEQ*NVAR], Type*const d2fds2_ptr[DIM*NEQ*NVAR*NVAR])
{
	compute_Flux_euler_0(flux_data,f_ptr);
	compute_Flux_euler_1(flux_data,dfds_ptr);
	compute_Flux_euler_2(flux_data,d2fds2_ptr);

	increment_pointers_T(DIM*NEQ,          (const Type**)f_ptr);
	increment_pointers_T(DIM*NEQ*NVAR,     (const Type**)dfds_ptr);
	increment_pointers_T(DIM*NEQ*NVAR*NVAR,(const Type**)d2fds2_ptr);
}

// Level 2 ********************************************************************************************************** //

static void compute_Flux_euler_0 (const struct Flux_Data_Euler*const flux_data, Type*const f_ptr[DIM*NEQ])
{
	const Type p = flux_data->p,
	           E = flux_data->E;

	IF_DIM_GE_1( const Type rhou = flux_data->rhouvw[0]; )
	IF_DIM_GE_2( const Type rhov = flux_data->rhouvw[1]; )
	IF_DIM_GE_3( const Type rhow = flux_data->rhouvw[2]; )

	IF_DIM_GE_1( const Type u = flux_data->uvw[0]; )
	IF_DIM_GE_2( const Type v = flux_data->uvw[1]; )
	IF_DIM_GE_3( const Type w = flux_data->uvw[2]; )

	int ind = 0;

	// f[:,0]
	IF_DIM_GE_1( *f_ptr[ind++] += rhou );
	IF_DIM_GE_2( *f_ptr[ind++] += rhov );
	IF_DIM_GE_3( *f_ptr[ind++] += rhow );

	// f[:,1]
	IF_DIM_GE_1( *f_ptr[ind++] += rhou*u + p );
	IF_DIM_GE_2( *f_ptr[ind++] += rhou*v );
	IF_DIM_GE_3( *f_ptr[ind++] += rhou*w );

	// f[:,2]
	IF_DIM_GE_2( *f_ptr[ind++] += rhov*u );
	IF_DIM_GE_2( *f_ptr[ind++] += rhov*v + p );
	IF_DIM_GE_3( *f_ptr[ind++] += rhov*w );

	// f[:,3]
	IF_DIM_GE_3( *f_ptr[ind++] += rhow*u );
	IF_DIM_GE_3( *f_ptr[ind++] += rhow*v );
	IF_DIM_GE_3( *f_ptr[ind++] += rhow*w + p );

	// f[:,4]
	IF_DIM_GE_1( *f_ptr[ind++] += (E+p)*u);
	IF_DIM_GE_2( *f_ptr[ind++] += (E+p)*v);
	IF_DIM_GE_3( *f_ptr[ind++] += (E+p)*w);
}

static void compute_Flux_euler_1 (const struct Flux_Data_Euler*const flux_data, Type*const dfds_ptr[DIM*NEQ*NVAR])
{
	const Type H     = flux_data->H,
	           alpha = flux_data->alpha,
	           beta  = flux_data->beta;

	IF_DIM_GE_1( const Type u = flux_data->uvw[0]; )
	IF_DIM_GE_2( const Type v = flux_data->uvw[1]; )
	IF_DIM_GE_3( const Type w = flux_data->uvw[2]; )

	int ind = 0;

	// dfds[:,:,0]
	// dfds[:,0,0]
	IF_DIM_GE_1( *dfds_ptr[ind++] += 0.0 );
	IF_DIM_GE_2( *dfds_ptr[ind++] += 0.0 );
	IF_DIM_GE_3( *dfds_ptr[ind++] += 0.0 );

	// dfds[:,1,0]
	IF_DIM_GE_1( *dfds_ptr[ind++] += -u*u + alpha );
	IF_DIM_GE_2( *dfds_ptr[ind++] += -u*v         );
	IF_DIM_GE_3( *dfds_ptr[ind++] += -u*w         );

	// dfds[:,2,0]
	IF_DIM_GE_2( *dfds_ptr[ind++] += -v*u         );
	IF_DIM_GE_2( *dfds_ptr[ind++] += -v*v + alpha );
	IF_DIM_GE_3( *dfds_ptr[ind++] += -v*w         );

	// dfds[:,3,0]
	IF_DIM_GE_3( *dfds_ptr[ind++] += -w*u         );
	IF_DIM_GE_3( *dfds_ptr[ind++] += -w*v         );
	IF_DIM_GE_3( *dfds_ptr[ind++] += -w*w + alpha );

	// dfds[:,4,0]
	IF_DIM_GE_1( *dfds_ptr[ind++] += u*beta );
	IF_DIM_GE_2( *dfds_ptr[ind++] += v*beta );
	IF_DIM_GE_3( *dfds_ptr[ind++] += w*beta );

	// dfds[:,:,1]
	// dfds[:,0,1]
	IF_DIM_GE_1( *dfds_ptr[ind++] += 1.0 );
	IF_DIM_GE_2( *dfds_ptr[ind++] += 0.0 );
	IF_DIM_GE_3( *dfds_ptr[ind++] += 0.0 );

	// dfds[:,1,1]
	IF_DIM_GE_1( *dfds_ptr[ind++] += -u*GM3 );
	IF_DIM_GE_2( *dfds_ptr[ind++] +=  v     );
	IF_DIM_GE_3( *dfds_ptr[ind++] +=  w     );

	// dfds[:,2,1]
	IF_DIM_GE_2( *dfds_ptr[ind++] +=  v     );
	IF_DIM_GE_2( *dfds_ptr[ind++] += -u*GM1 );
	IF_DIM_GE_3( *dfds_ptr[ind++] +=  0.0   );

	// dfds[:,3,1]
	IF_DIM_GE_3( *dfds_ptr[ind++] +=  w     );
	IF_DIM_GE_3( *dfds_ptr[ind++] +=  0.0   );
	IF_DIM_GE_3( *dfds_ptr[ind++] += -u*GM1 );

	// dfds[:,4,1]
	IF_DIM_GE_1( *dfds_ptr[ind++] += -GM1*u*u + H );
	IF_DIM_GE_2( *dfds_ptr[ind++] += -GM1*u*v     );
	IF_DIM_GE_3( *dfds_ptr[ind++] += -GM1*u*w     );

	// dfds[:,:,2]
	// dfds[:,0,2]
	IF_DIM_GE_2( *dfds_ptr[ind++] += 0.0 );
	IF_DIM_GE_2( *dfds_ptr[ind++] += 1.0 );
	IF_DIM_GE_3( *dfds_ptr[ind++] += 0.0 );

	// dfds[:,1,2]
	IF_DIM_GE_2( *dfds_ptr[ind++] += -v*GM1 );
	IF_DIM_GE_2( *dfds_ptr[ind++] +=  u     );
	IF_DIM_GE_3( *dfds_ptr[ind++] +=  0.0   );

	// dfds[:,2,2]
	IF_DIM_GE_2( *dfds_ptr[ind++] +=  u     );
	IF_DIM_GE_2( *dfds_ptr[ind++] += -v*GM3 );
	IF_DIM_GE_3( *dfds_ptr[ind++] +=  w     );

	// dfds[:,3,2]
	IF_DIM_GE_3( *dfds_ptr[ind++] +=  0.0   );
	IF_DIM_GE_3( *dfds_ptr[ind++] +=  w     );
	IF_DIM_GE_3( *dfds_ptr[ind++] += -v*GM1 );

	// dfds[:,4,2]
	IF_DIM_GE_2( *dfds_ptr[ind++] += -GM1*v*u     );
	IF_DIM_GE_2( *dfds_ptr[ind++] += -GM1*v*v + H );
	IF_DIM_GE_3( *dfds_ptr[ind++] += -GM1*v*w     );

	// dfds[:,:,3]
	// dfds[:,0,3]
	IF_DIM_GE_3( *dfds_ptr[ind++] += 0.0 );
	IF_DIM_GE_3( *dfds_ptr[ind++] += 0.0 );
	IF_DIM_GE_3( *dfds_ptr[ind++] += 1.0 );

	// dfds[:,1,3]
	IF_DIM_GE_3( *dfds_ptr[ind++] += -w*GM1 );
	IF_DIM_GE_3( *dfds_ptr[ind++] +=  0.0   );
	IF_DIM_GE_3( *dfds_ptr[ind++] +=  u     );

	// dfds[:,2,3]
	IF_DIM_GE_3( *dfds_ptr[ind++] +=  0.0   );
	IF_DIM_GE_3( *dfds_ptr[ind++] += -w*GM1 );
	IF_DIM_GE_3( *dfds_ptr[ind++] +=  v     );

	// dfds[:,3,3]
	IF_DIM_GE_3( *dfds_ptr[ind++] +=  u     );
	IF_DIM_GE_3( *dfds_ptr[ind++] +=  v     );
	IF_DIM_GE_3( *dfds_ptr[ind++] += -w*GM3 );

	// dfds[:,4,3]
	IF_DIM_GE_3( *dfds_ptr[ind++] += -GM1*w*u     );
	IF_DIM_GE_3( *dfds_ptr[ind++] += -GM1*w*v     );
	IF_DIM_GE_3( *dfds_ptr[ind++] += -GM1*w*w + H );

	// dfds[:,:,4]
	// dfds[:,0,4]
	IF_DIM_GE_1( *dfds_ptr[ind++] += 0.0 );
	IF_DIM_GE_2( *dfds_ptr[ind++] += 0.0 );
	IF_DIM_GE_3( *dfds_ptr[ind++] += 0.0 );

	// dfds[:,1,4]
	IF_DIM_GE_1( *dfds_ptr[ind++] += GM1 );
	IF_DIM_GE_2( *dfds_ptr[ind++] += 0.0 );
	IF_DIM_GE_3( *dfds_ptr[ind++] += 0.0 );

	// dfds[:,2,4]
	IF_DIM_GE_2( *dfds_ptr[ind++] += 0.0 );
	IF_DIM_GE_2( *dfds_ptr[ind++] += GM1 );
	IF_DIM_GE_3( *dfds_ptr[ind++] += 0.0 );

	// dfds[:,3,4]
	IF_DIM_GE_3( *dfds_ptr[ind++] += 0.0 );
	IF_DIM_GE_3( *dfds_ptr[ind++] += 0.0 );
	IF_DIM_GE_3( *dfds_ptr[ind++] += GM1 );

	// dfds[:,4,4]
	IF_DIM_GE_1( *dfds_ptr[ind++] += u*GAMMA );
	IF_DIM_GE_2( *dfds_ptr[ind++] += v*GAMMA );
	IF_DIM_GE_3( *dfds_ptr[ind++] += w*GAMMA );
}

static void compute_Flux_euler_2
	(const struct Flux_Data_Euler*const flux_data, Type*const d2fds2_ptr[DIM*NEQ*NVAR*NVAR])
{
	const Type rho_inv = flux_data->rho_inv,
	           E       = flux_data->E,
	           V2      = flux_data->V2,
	           alpha   = flux_data->alpha,
	           beta    = flux_data->beta;

	IF_DIM_GE_1( const Type u = flux_data->uvw[0]; )
	IF_DIM_GE_2( const Type v = flux_data->uvw[1]; )
	IF_DIM_GE_3( const Type w = flux_data->uvw[2]; )

	const Type rho2_inv = rho_inv*rho_inv;

	Type dalpha_dW[NVAR];
	              dalpha_dW[0]     = -2.0*rho_inv*alpha;
	IF_DIM_GE_1(  dalpha_dW[1]     =  GM1*rho_inv*u );
	IF_DIM_GE_2(  dalpha_dW[2]     =  GM1*rho_inv*v );
	IF_DIM_GE_3(  dalpha_dW[3]     =  GM1*rho_inv*w );
	              dalpha_dW[DIM+1] =  0.0;

	Type dH_dW[NVAR];
	             dH_dW[0]     = -GAMMA*rho2_inv*E + GM1*rho_inv*V2;
	IF_DIM_GE_1( dH_dW[1]     = -GM1*rho_inv*u );
	IF_DIM_GE_2( dH_dW[2]     = -GM1*rho_inv*v );
	IF_DIM_GE_3( dH_dW[3]     = -GM1*rho_inv*w );
	             dH_dW[DIM+1] =  GAMMA*rho_inv;

	Type dbeta_dW[NVAR];
	for (int i = 0; i < NVAR; ++i)
		dbeta_dW[i] = dalpha_dW[i]-dH_dW[i];

	int ind = 0;

	// d2fds2[:,:,:,0]
	// d2fds2[:,:,0,0]
	// d2fds2[:,0,0,0]
	IF_DIM_GE_1( *d2fds2_ptr[ind++] += 0.0 );
	IF_DIM_GE_2( *d2fds2_ptr[ind++] += 0.0 );
	IF_DIM_GE_3( *d2fds2_ptr[ind++] += 0.0 );

	// d2fds2[:,1,0,0]
	IF_DIM_GE_1( *d2fds2_ptr[ind++] += rho_inv*(2.0*u*u - GM1*V2) );
	IF_DIM_GE_2( *d2fds2_ptr[ind++] += rho_inv*(2.0*u*v         ) );
	IF_DIM_GE_3( *d2fds2_ptr[ind++] += rho_inv*(2.0*u*w         ) );

	// d2fds2[:,2,0,0]
	IF_DIM_GE_2( *d2fds2_ptr[ind++] += rho_inv*(2.0*v*u         ) );
	IF_DIM_GE_2( *d2fds2_ptr[ind++] += rho_inv*(2.0*v*v - GM1*V2) );
	IF_DIM_GE_3( *d2fds2_ptr[ind++] += rho_inv*(2.0*v*w         ) );

	// d2fds2[:,3,0,0]
	IF_DIM_GE_3( *d2fds2_ptr[ind++] += rho_inv*(2.0*w*u         ) );
	IF_DIM_GE_3( *d2fds2_ptr[ind++] += rho_inv*(2.0*w*v         ) );
	IF_DIM_GE_3( *d2fds2_ptr[ind++] += rho_inv*(2.0*w*w - GM1*V2) );

	// d2fds2[:,4,0,0]
	IF_DIM_GE_1( *d2fds2_ptr[ind++] += u*dbeta_dW[0] - rho_inv*u*beta );
	IF_DIM_GE_2( *d2fds2_ptr[ind++] += v*dbeta_dW[0] - rho_inv*v*beta );
	IF_DIM_GE_3( *d2fds2_ptr[ind++] += w*dbeta_dW[0] - rho_inv*w*beta );

	// d2fds2[:,:,1,0]
	// d2fds2[:,0,1,0]
	IF_DIM_GE_1( *d2fds2_ptr[ind++] += 0.0 );
	IF_DIM_GE_2( *d2fds2_ptr[ind++] += 0.0 );
	IF_DIM_GE_3( *d2fds2_ptr[ind++] += 0.0 );

	// d2fds2[:,1,1,0]
	IF_DIM_GE_1( *d2fds2_ptr[ind++] +=  rho_inv*u*GM3 );
	IF_DIM_GE_2( *d2fds2_ptr[ind++] += -rho_inv*v     );
	IF_DIM_GE_3( *d2fds2_ptr[ind++] += -rho_inv*w     );

	// d2fds2[:,2,1,0]
	IF_DIM_GE_2( *d2fds2_ptr[ind++] += -rho_inv*v     );
	IF_DIM_GE_2( *d2fds2_ptr[ind++] +=  rho_inv*u*GM1 );
	IF_DIM_GE_3( *d2fds2_ptr[ind++] +=  0.0           );

	// d2fds2[:,3,1,0]
	IF_DIM_GE_3( *d2fds2_ptr[ind++] += -rho_inv*w     );
	IF_DIM_GE_3( *d2fds2_ptr[ind++] +=  0.0           );
	IF_DIM_GE_3( *d2fds2_ptr[ind++] +=  rho_inv*u*GM1 );

	// d2fds2[:,4,1,0]
	IF_DIM_GE_1( *d2fds2_ptr[ind++] += 2.0*GM1*rho_inv*u*u + dH_dW[0] );
	IF_DIM_GE_2( *d2fds2_ptr[ind++] += 2.0*GM1*rho_inv*u*v            );
	IF_DIM_GE_3( *d2fds2_ptr[ind++] += 2.0*GM1*rho_inv*u*w            );

	// d2fds2[:,:,2,0]
	// d2fds2[:,0,2,0]
	IF_DIM_GE_1( *d2fds2_ptr[ind++] += 0.0 );
	IF_DIM_GE_2( *d2fds2_ptr[ind++] += 0.0 );
	IF_DIM_GE_3( *d2fds2_ptr[ind++] += 0.0 );

	// d2fds2[:,1,2,0]
	IF_DIM_GE_2( *d2fds2_ptr[ind++] +=  rho_inv*v*GM1 );
	IF_DIM_GE_2( *d2fds2_ptr[ind++] += -rho_inv*u     );
	IF_DIM_GE_3( *d2fds2_ptr[ind++] +=  0.0           );

	// d2fds2[:,2,2,0]
	IF_DIM_GE_2( *d2fds2_ptr[ind++] += -rho_inv*u     );
	IF_DIM_GE_2( *d2fds2_ptr[ind++] +=  rho_inv*v*GM3 );
	IF_DIM_GE_3( *d2fds2_ptr[ind++] += -rho_inv*w     );

	// d2fds2[:,3,2,0]
	IF_DIM_GE_3( *d2fds2_ptr[ind++] +=  0.0           );
	IF_DIM_GE_3( *d2fds2_ptr[ind++] += -rho_inv*w     );
	IF_DIM_GE_3( *d2fds2_ptr[ind++] +=  rho_inv*v*GM1 );

	// d2fds2[:,4,2,0]
	IF_DIM_GE_2( *d2fds2_ptr[ind++] += 2.0*GM1*rho_inv*v*u            );
	IF_DIM_GE_2( *d2fds2_ptr[ind++] += 2.0*GM1*rho_inv*v*v + dH_dW[0] );
	IF_DIM_GE_3( *d2fds2_ptr[ind++] += 2.0*GM1*rho_inv*v*w            );

	// d2fds2[:,:,3,0]
	// d2fds2[:,0,3,0]
	IF_DIM_GE_3( *d2fds2_ptr[ind++] += 0.0 );
	IF_DIM_GE_3( *d2fds2_ptr[ind++] += 0.0 );
	IF_DIM_GE_3( *d2fds2_ptr[ind++] += 0.0 );

	// d2fds2[:,1,3,0]
	IF_DIM_GE_3( *d2fds2_ptr[ind++] +=  rho_inv*w*GM1 );
	IF_DIM_GE_3( *d2fds2_ptr[ind++] +=  0.0           );
	IF_DIM_GE_3( *d2fds2_ptr[ind++] += -rho_inv*u     );

	// d2fds2[:,2,3,0]
	IF_DIM_GE_3( *d2fds2_ptr[ind++] +=  0.0           );
	IF_DIM_GE_3( *d2fds2_ptr[ind++] +=  rho_inv*w*GM1 );
	IF_DIM_GE_3( *d2fds2_ptr[ind++] += -rho_inv*v     );

	// d2fds2[:,3,3,0]
	IF_DIM_GE_3( *d2fds2_ptr[ind++] += -rho_inv*u     );
	IF_DIM_GE_3( *d2fds2_ptr[ind++] += -rho_inv*v     );
	IF_DIM_GE_3( *d2fds2_ptr[ind++] +=  rho_inv*w*GM3 );

	// d2fds2[:,4,3,0]
	IF_DIM_GE_3( *d2fds2_ptr[ind++] += 2.0*GM1*rho_inv*w*u            );
	IF_DIM_GE_3( *d2fds2_ptr[ind++] += 2.0*GM1*rho_inv*w*v            );
	IF_DIM_GE_3( *d2fds2_ptr[ind++] += 2.0*GM1*rho_inv*w*w + dH_dW[0] );

	// d2fds2[:,:,4,0]
	// d2fds2[:,0,4,0]
	IF_DIM_GE_1( *d2fds2_ptr[ind++] += 0.0 );
	IF_DIM_GE_2( *d2fds2_ptr[ind++] += 0.0 );
	IF_DIM_GE_3( *d2fds2_ptr[ind++] += 0.0 );

	// d2fds2[:,1,4,0]
	IF_DIM_GE_1( *d2fds2_ptr[ind++] += 0.0 );
	IF_DIM_GE_2( *d2fds2_ptr[ind++] += 0.0 );
	IF_DIM_GE_3( *d2fds2_ptr[ind++] += 0.0 );

	// d2fds2[:,2,4,0]
	IF_DIM_GE_2( *d2fds2_ptr[ind++] += 0.0 );
	IF_DIM_GE_2( *d2fds2_ptr[ind++] += 0.0 );
	IF_DIM_GE_3( *d2fds2_ptr[ind++] += 0.0 );

	// d2fds2[:,3,4,0]
	IF_DIM_GE_3( *d2fds2_ptr[ind++] += 0.0 );
	IF_DIM_GE_3( *d2fds2_ptr[ind++] += 0.0 );
	IF_DIM_GE_3( *d2fds2_ptr[ind++] += 0.0 );

	// d2fds2[:,4,4,0]
	IF_DIM_GE_1( *d2fds2_ptr[ind++] += -GAMMA*rho_inv*u );
	IF_DIM_GE_2( *d2fds2_ptr[ind++] += -GAMMA*rho_inv*v );
	IF_DIM_GE_3( *d2fds2_ptr[ind++] += -GAMMA*rho_inv*w );

	// d2fds2[:,:,:,1]
	// d2fds2[:,:,0,1]
	// d2fds2[:,0,0,1]
	IF_DIM_GE_1( *d2fds2_ptr[ind++] += 0.0 );
	IF_DIM_GE_2( *d2fds2_ptr[ind++] += 0.0 );
	IF_DIM_GE_3( *d2fds2_ptr[ind++] += 0.0 );

	// d2fds2[:,1,0,1]
	IF_DIM_GE_1( *d2fds2_ptr[ind++] +=  u*rho_inv*GM3 );
	IF_DIM_GE_2( *d2fds2_ptr[ind++] += -v*rho_inv     );
	IF_DIM_GE_3( *d2fds2_ptr[ind++] += -w*rho_inv     );

	// d2fds2[:,2,0,1]
	IF_DIM_GE_2( *d2fds2_ptr[ind++] += -v*rho_inv     );
	IF_DIM_GE_2( *d2fds2_ptr[ind++] +=  u*rho_inv*GM1 );
	IF_DIM_GE_3( *d2fds2_ptr[ind++] +=  0.0           );

	// d2fds2[:,3,0,1]
	IF_DIM_GE_3( *d2fds2_ptr[ind++] += -w*rho_inv     );
	IF_DIM_GE_3( *d2fds2_ptr[ind++] +=  0.0           );
	IF_DIM_GE_3( *d2fds2_ptr[ind++] +=  u*rho_inv*GM1 );

	// d2fds2[:,4,0,1]
	IF_DIM_GE_1( *d2fds2_ptr[ind++] += u*dbeta_dW[1] + rho_inv*beta );
	IF_DIM_GE_2( *d2fds2_ptr[ind++] += v*dbeta_dW[1]                );
	IF_DIM_GE_3( *d2fds2_ptr[ind++] += w*dbeta_dW[1]                );

	// d2fds2[:,:,1,1]
	// d2fds2[:,0,1,1]
	IF_DIM_GE_1( *d2fds2_ptr[ind++] += 0.0 );
	IF_DIM_GE_2( *d2fds2_ptr[ind++] += 0.0 );
	IF_DIM_GE_3( *d2fds2_ptr[ind++] += 0.0 );

	// d2fds2[:,1,1,1]
	IF_DIM_GE_1( *d2fds2_ptr[ind++] += -rho_inv*GM3 );
	IF_DIM_GE_2( *d2fds2_ptr[ind++] +=  0.0         );
	IF_DIM_GE_3( *d2fds2_ptr[ind++] +=  0.0         );

	// d2fds2[:,2,1,1]
	IF_DIM_GE_2( *d2fds2_ptr[ind++] +=  0.0         );
	IF_DIM_GE_2( *d2fds2_ptr[ind++] += -rho_inv*GM1 );
	IF_DIM_GE_3( *d2fds2_ptr[ind++] +=  0.0         );

	// d2fds2[:,3,1,1]
	IF_DIM_GE_3( *d2fds2_ptr[ind++] +=  0.0         );
	IF_DIM_GE_3( *d2fds2_ptr[ind++] +=  0.0         );
	IF_DIM_GE_3( *d2fds2_ptr[ind++] += -rho_inv*GM1 );

	// d2fds2[:,4,1,1]
	IF_DIM_GE_1( *d2fds2_ptr[ind++] += -GM1*rho_inv*u*3.0 );
	IF_DIM_GE_2( *d2fds2_ptr[ind++] += -GM1*rho_inv*v     );
	IF_DIM_GE_3( *d2fds2_ptr[ind++] += -GM1*rho_inv*w     );

	// d2fds2[:,:,2,1]
	// d2fds2[:,0,2,1]
	IF_DIM_GE_1( *d2fds2_ptr[ind++] += 0.0 );
	IF_DIM_GE_2( *d2fds2_ptr[ind++] += 0.0 );
	IF_DIM_GE_3( *d2fds2_ptr[ind++] += 0.0 );

	// d2fds2[:,1,2,1]
	IF_DIM_GE_2( *d2fds2_ptr[ind++] +=  0.0         );
	IF_DIM_GE_2( *d2fds2_ptr[ind++] +=  rho_inv     );
	IF_DIM_GE_3( *d2fds2_ptr[ind++] +=  0.0         );

	// d2fds2[:,2,2,1]
	IF_DIM_GE_2( *d2fds2_ptr[ind++] +=  rho_inv     );
	IF_DIM_GE_2( *d2fds2_ptr[ind++] +=  0.0         );
	IF_DIM_GE_3( *d2fds2_ptr[ind++] +=  0.0         );

	// d2fds2[:,3,2,1]
	IF_DIM_GE_3( *d2fds2_ptr[ind++] +=  0.0         );
	IF_DIM_GE_3( *d2fds2_ptr[ind++] +=  0.0         );
	IF_DIM_GE_3( *d2fds2_ptr[ind++] +=  0.0         );

	// d2fds2[:,4,2,1]
	IF_DIM_GE_2( *d2fds2_ptr[ind++] += -GM1*rho_inv*v     );
	IF_DIM_GE_2( *d2fds2_ptr[ind++] += -GM1*rho_inv*u     );
	IF_DIM_GE_3( *d2fds2_ptr[ind++] +=  0.0               );

	// d2fds2[:,:,3,1]
	// d2fds2[:,0,3,1]
	IF_DIM_GE_3( *d2fds2_ptr[ind++] += 0.0 );
	IF_DIM_GE_3( *d2fds2_ptr[ind++] += 0.0 );
	IF_DIM_GE_3( *d2fds2_ptr[ind++] += 0.0 );

	// d2fds2[:,1,3,1]
	IF_DIM_GE_3( *d2fds2_ptr[ind++] +=  0.0         );
	IF_DIM_GE_3( *d2fds2_ptr[ind++] +=  0.0         );
	IF_DIM_GE_3( *d2fds2_ptr[ind++] +=  rho_inv     );

	// d2fds2[:,2,3,1]
	IF_DIM_GE_3( *d2fds2_ptr[ind++] +=  0.0         );
	IF_DIM_GE_3( *d2fds2_ptr[ind++] +=  0.0         );
	IF_DIM_GE_3( *d2fds2_ptr[ind++] +=  0.0         );

	// d2fds2[:,3,3,1]
	IF_DIM_GE_3( *d2fds2_ptr[ind++] +=  rho_inv     );
	IF_DIM_GE_3( *d2fds2_ptr[ind++] +=  0.0         );
	IF_DIM_GE_3( *d2fds2_ptr[ind++] +=  0.0         );

	// d2fds2[:,4,3,1]
	IF_DIM_GE_3( *d2fds2_ptr[ind++] += -GM1*rho_inv*w     );
	IF_DIM_GE_3( *d2fds2_ptr[ind++] +=  0.0               );
	IF_DIM_GE_3( *d2fds2_ptr[ind++] += -GM1*rho_inv*u     );

	// d2fds2[:,:,4,1]
	// d2fds2[:,0,4,1]
	IF_DIM_GE_1( *d2fds2_ptr[ind++] += 0.0 );
	IF_DIM_GE_2( *d2fds2_ptr[ind++] += 0.0 );
	IF_DIM_GE_3( *d2fds2_ptr[ind++] += 0.0 );

	// d2fds2[:,1,4,1]
	IF_DIM_GE_1( *d2fds2_ptr[ind++] += 0.0 );
	IF_DIM_GE_2( *d2fds2_ptr[ind++] += 0.0 );
	IF_DIM_GE_3( *d2fds2_ptr[ind++] += 0.0 );

	// d2fds2[:,2,4,1]
	IF_DIM_GE_2( *d2fds2_ptr[ind++] += 0.0 );
	IF_DIM_GE_2( *d2fds2_ptr[ind++] += 0.0 );
	IF_DIM_GE_3( *d2fds2_ptr[ind++] += 0.0 );

	// d2fds2[:,3,4,1]
	IF_DIM_GE_3( *d2fds2_ptr[ind++] += 0.0 );
	IF_DIM_GE_3( *d2fds2_ptr[ind++] += 0.0 );
	IF_DIM_GE_3( *d2fds2_ptr[ind++] += 0.0 );

	// d2fds2[:,4,4,1]
	IF_DIM_GE_1( *d2fds2_ptr[ind++] += GAMMA*rho_inv );
	IF_DIM_GE_2( *d2fds2_ptr[ind++] += 0.0           );
	IF_DIM_GE_3( *d2fds2_ptr[ind++] += 0.0           );

	// d2fds2[:,:,:,2]
	// d2fds2[:,:,0,2]
	// d2fds2[:,0,0,2]
	IF_DIM_GE_2( *d2fds2_ptr[ind++] += 0.0 );
	IF_DIM_GE_2( *d2fds2_ptr[ind++] += 0.0 );
	IF_DIM_GE_3( *d2fds2_ptr[ind++] += 0.0 );

	// d2fds2[:,1,0,2]
	IF_DIM_GE_2( *d2fds2_ptr[ind++] +=  v*rho_inv*GM1 );
	IF_DIM_GE_2( *d2fds2_ptr[ind++] += -u*rho_inv     );
	IF_DIM_GE_3( *d2fds2_ptr[ind++] +=  0.0           );

	// d2fds2[:,2,0,2]
	IF_DIM_GE_2( *d2fds2_ptr[ind++] += -u*rho_inv     );
	IF_DIM_GE_2( *d2fds2_ptr[ind++] +=  v*rho_inv*GM3 );
	IF_DIM_GE_3( *d2fds2_ptr[ind++] += -w*rho_inv     );

	// d2fds2[:,3,0,2]
	IF_DIM_GE_3( *d2fds2_ptr[ind++] +=  0.0           );
	IF_DIM_GE_3( *d2fds2_ptr[ind++] += -w*rho_inv     );
	IF_DIM_GE_3( *d2fds2_ptr[ind++] +=  v*rho_inv*GM1 );

	// d2fds2[:,4,0,2]
	IF_DIM_GE_2( *d2fds2_ptr[ind++] += u*dbeta_dW[2]                );
	IF_DIM_GE_2( *d2fds2_ptr[ind++] += v*dbeta_dW[2] + rho_inv*beta );
	IF_DIM_GE_3( *d2fds2_ptr[ind++] += w*dbeta_dW[2]                );

	// d2fds2[:,:,1,2]
	// d2fds2[:,0,1,2]
	IF_DIM_GE_2( *d2fds2_ptr[ind++] += 0.0 );
	IF_DIM_GE_2( *d2fds2_ptr[ind++] += 0.0 );
	IF_DIM_GE_3( *d2fds2_ptr[ind++] += 0.0 );

	// d2fds2[:,1,1,2]
	IF_DIM_GE_2( *d2fds2_ptr[ind++] +=  0.0         );
	IF_DIM_GE_2( *d2fds2_ptr[ind++] +=  rho_inv     );
	IF_DIM_GE_3( *d2fds2_ptr[ind++] +=  0.0         );

	// d2fds2[:,2,1,2]
	IF_DIM_GE_2( *d2fds2_ptr[ind++] +=  rho_inv     );
	IF_DIM_GE_2( *d2fds2_ptr[ind++] +=  0.0         );
	IF_DIM_GE_3( *d2fds2_ptr[ind++] +=  0.0         );

	// d2fds2[:,3,1,2]
	IF_DIM_GE_3( *d2fds2_ptr[ind++] +=  0.0         );
	IF_DIM_GE_3( *d2fds2_ptr[ind++] +=  0.0         );
	IF_DIM_GE_3( *d2fds2_ptr[ind++] +=  0.0         );

	// d2fds2[:,4,1,2]
	IF_DIM_GE_2( *d2fds2_ptr[ind++] += -GM1*rho_inv*v     );
	IF_DIM_GE_2( *d2fds2_ptr[ind++] += -GM1*rho_inv*u     );
	IF_DIM_GE_3( *d2fds2_ptr[ind++] +=  0.0               );

	// d2fds2[:,:,2,2]
	// d2fds2[:,0,2,2]
	IF_DIM_GE_2( *d2fds2_ptr[ind++] += 0.0 );
	IF_DIM_GE_2( *d2fds2_ptr[ind++] += 0.0 );
	IF_DIM_GE_3( *d2fds2_ptr[ind++] += 0.0 );

	// d2fds2[:,1,2,2]
	IF_DIM_GE_2( *d2fds2_ptr[ind++] += -rho_inv*GM1 );
	IF_DIM_GE_2( *d2fds2_ptr[ind++] +=  0.0         );
	IF_DIM_GE_3( *d2fds2_ptr[ind++] +=  0.0         );

	// d2fds2[:,2,2,2]
	IF_DIM_GE_2( *d2fds2_ptr[ind++] +=  0.0         );
	IF_DIM_GE_2( *d2fds2_ptr[ind++] += -rho_inv*GM3 );
	IF_DIM_GE_3( *d2fds2_ptr[ind++] +=  0.0         );

	// d2fds2[:,3,2,2]
	IF_DIM_GE_3( *d2fds2_ptr[ind++] +=  0.0         );
	IF_DIM_GE_3( *d2fds2_ptr[ind++] +=  0.0         );
	IF_DIM_GE_3( *d2fds2_ptr[ind++] += -rho_inv*GM1 );

	// d2fds2[:,4,2,2]
	IF_DIM_GE_2( *d2fds2_ptr[ind++] += -GM1*rho_inv*u     );
	IF_DIM_GE_2( *d2fds2_ptr[ind++] += -GM1*rho_inv*v*3.0 );
	IF_DIM_GE_3( *d2fds2_ptr[ind++] += -GM1*rho_inv*w     );

	// d2fds2[:,:,3,2]
	// d2fds2[:,0,3,2]
	IF_DIM_GE_3( *d2fds2_ptr[ind++] += 0.0 );
	IF_DIM_GE_3( *d2fds2_ptr[ind++] += 0.0 );
	IF_DIM_GE_3( *d2fds2_ptr[ind++] += 0.0 );

	// d2fds2[:,1,3,2]
	IF_DIM_GE_3( *d2fds2_ptr[ind++] +=  0.0         );
	IF_DIM_GE_3( *d2fds2_ptr[ind++] +=  0.0         );
	IF_DIM_GE_3( *d2fds2_ptr[ind++] +=  0.0         );

	// d2fds2[:,2,3,2]
	IF_DIM_GE_3( *d2fds2_ptr[ind++] +=  0.0         );
	IF_DIM_GE_3( *d2fds2_ptr[ind++] +=  0.0         );
	IF_DIM_GE_3( *d2fds2_ptr[ind++] +=  rho_inv     );

	// d2fds2[:,3,3,2]
	IF_DIM_GE_3( *d2fds2_ptr[ind++] +=  0.0         );
	IF_DIM_GE_3( *d2fds2_ptr[ind++] +=  rho_inv     );
	IF_DIM_GE_3( *d2fds2_ptr[ind++] +=  0.0         );

	// d2fds2[:,4,3,2]
	IF_DIM_GE_3( *d2fds2_ptr[ind++] +=  0.0               );
	IF_DIM_GE_3( *d2fds2_ptr[ind++] += -GM1*rho_inv*v     );
	IF_DIM_GE_3( *d2fds2_ptr[ind++] += -GM1*rho_inv*v     );

	// d2fds2[:,:,4,2]
	// d2fds2[:,0,4,2]
	IF_DIM_GE_2( *d2fds2_ptr[ind++] += 0.0 );
	IF_DIM_GE_2( *d2fds2_ptr[ind++] += 0.0 );
	IF_DIM_GE_3( *d2fds2_ptr[ind++] += 0.0 );

	// d2fds2[:,1,4,2]
	IF_DIM_GE_2( *d2fds2_ptr[ind++] += 0.0 );
	IF_DIM_GE_2( *d2fds2_ptr[ind++] += 0.0 );
	IF_DIM_GE_3( *d2fds2_ptr[ind++] += 0.0 );

	// d2fds2[:,2,4,2]
	IF_DIM_GE_2( *d2fds2_ptr[ind++] += 0.0 );
	IF_DIM_GE_2( *d2fds2_ptr[ind++] += 0.0 );
	IF_DIM_GE_3( *d2fds2_ptr[ind++] += 0.0 );

	// d2fds2[:,3,4,2]
	IF_DIM_GE_3( *d2fds2_ptr[ind++] += 0.0 );
	IF_DIM_GE_3( *d2fds2_ptr[ind++] += 0.0 );
	IF_DIM_GE_3( *d2fds2_ptr[ind++] += 0.0 );

	// d2fds2[:,4,4,2]
	IF_DIM_GE_2( *d2fds2_ptr[ind++] += 0.0           );
	IF_DIM_GE_2( *d2fds2_ptr[ind++] += GAMMA*rho_inv );
	IF_DIM_GE_3( *d2fds2_ptr[ind++] += 0.0           );

	// d2fds2[:,:,:,3]
	// d2fds2[:,:,0,3]
	// d2fds2[:,0,0,3]
	IF_DIM_GE_3( *d2fds2_ptr[ind++] += 0.0 );
	IF_DIM_GE_3( *d2fds2_ptr[ind++] += 0.0 );
	IF_DIM_GE_3( *d2fds2_ptr[ind++] += 0.0 );

	// d2fds2[:,1,0,3]
	IF_DIM_GE_3( *d2fds2_ptr[ind++] +=  w*rho_inv*GM1 );
	IF_DIM_GE_3( *d2fds2_ptr[ind++] +=  0.0           );
	IF_DIM_GE_3( *d2fds2_ptr[ind++] += -u*rho_inv     );

	// d2fds2[:,2,0,3]
	IF_DIM_GE_3( *d2fds2_ptr[ind++] +=  0.0           );
	IF_DIM_GE_3( *d2fds2_ptr[ind++] +=  w*rho_inv*GM1 );
	IF_DIM_GE_3( *d2fds2_ptr[ind++] += -v*rho_inv     );

	// d2fds2[:,3,0,3]
	IF_DIM_GE_3( *d2fds2_ptr[ind++] += -u*rho_inv     );
	IF_DIM_GE_3( *d2fds2_ptr[ind++] += -v*rho_inv     );
	IF_DIM_GE_3( *d2fds2_ptr[ind++] +=  w*rho_inv*GM3 );

	// d2fds2[:,4,0,3]
	IF_DIM_GE_3( *d2fds2_ptr[ind++] += u*dbeta_dW[3]                );
	IF_DIM_GE_3( *d2fds2_ptr[ind++] += v*dbeta_dW[3]                );
	IF_DIM_GE_3( *d2fds2_ptr[ind++] += w*dbeta_dW[3] + rho_inv*beta );

	// d2fds2[:,:,1,3]
	// d2fds2[:,0,1,3]
	IF_DIM_GE_3( *d2fds2_ptr[ind++] += 0.0 );
	IF_DIM_GE_3( *d2fds2_ptr[ind++] += 0.0 );
	IF_DIM_GE_3( *d2fds2_ptr[ind++] += 0.0 );

	// d2fds2[:,1,1,3]
	IF_DIM_GE_3( *d2fds2_ptr[ind++] +=  0.0         );
	IF_DIM_GE_3( *d2fds2_ptr[ind++] +=  0.0         );
	IF_DIM_GE_3( *d2fds2_ptr[ind++] +=  rho_inv     );

	// d2fds2[:,2,1,3]
	IF_DIM_GE_3( *d2fds2_ptr[ind++] +=  0.0         );
	IF_DIM_GE_3( *d2fds2_ptr[ind++] +=  0.0         );
	IF_DIM_GE_3( *d2fds2_ptr[ind++] +=  0.0         );

	// d2fds2[:,3,1,3]
	IF_DIM_GE_3( *d2fds2_ptr[ind++] +=  rho_inv     );
	IF_DIM_GE_3( *d2fds2_ptr[ind++] +=  0.0         );
	IF_DIM_GE_3( *d2fds2_ptr[ind++] +=  0.0         );

	// d2fds2[:,4,1,3]
	IF_DIM_GE_3( *d2fds2_ptr[ind++] += -GM1*rho_inv*v     );
	IF_DIM_GE_3( *d2fds2_ptr[ind++] +=  0.0               );
	IF_DIM_GE_3( *d2fds2_ptr[ind++] += -GM1*rho_inv*u     );

	// d2fds2[:,:,2,3]
	// d2fds2[:,0,2,3]
	IF_DIM_GE_3( *d2fds2_ptr[ind++] += 0.0 );
	IF_DIM_GE_3( *d2fds2_ptr[ind++] += 0.0 );
	IF_DIM_GE_3( *d2fds2_ptr[ind++] += 0.0 );

	// d2fds2[:,1,2,3]
	IF_DIM_GE_3( *d2fds2_ptr[ind++] +=  0.0         );
	IF_DIM_GE_3( *d2fds2_ptr[ind++] +=  0.0         );
	IF_DIM_GE_3( *d2fds2_ptr[ind++] +=  0.0         );

	// d2fds2[:,2,2,3]
	IF_DIM_GE_3( *d2fds2_ptr[ind++] +=  0.0         );
	IF_DIM_GE_3( *d2fds2_ptr[ind++] +=  0.0         );
	IF_DIM_GE_3( *d2fds2_ptr[ind++] +=  rho_inv     );

	// d2fds2[:,3,2,3]
	IF_DIM_GE_3( *d2fds2_ptr[ind++] +=  0.0         );
	IF_DIM_GE_3( *d2fds2_ptr[ind++] +=  rho_inv     );
	IF_DIM_GE_3( *d2fds2_ptr[ind++] +=  0.0         );

	// d2fds2[:,4,2,3]
	IF_DIM_GE_3( *d2fds2_ptr[ind++] +=  0.0               );
	IF_DIM_GE_3( *d2fds2_ptr[ind++] += -GM1*rho_inv*w     );
	IF_DIM_GE_3( *d2fds2_ptr[ind++] += -GM1*rho_inv*v     );

	// d2fds2[:,:,3,3]
	// d2fds2[:,0,3,3]
	IF_DIM_GE_3( *d2fds2_ptr[ind++] += 0.0 );
	IF_DIM_GE_3( *d2fds2_ptr[ind++] += 0.0 );
	IF_DIM_GE_3( *d2fds2_ptr[ind++] += 0.0 );

	// d2fds2[:,1,3,3]
	IF_DIM_GE_3( *d2fds2_ptr[ind++] += -rho_inv*GM1 );
	IF_DIM_GE_3( *d2fds2_ptr[ind++] +=  0.0         );
	IF_DIM_GE_3( *d2fds2_ptr[ind++] +=  0.0         );

	// d2fds2[:,2,3,3]
	IF_DIM_GE_3( *d2fds2_ptr[ind++] +=  0.0         );
	IF_DIM_GE_3( *d2fds2_ptr[ind++] += -rho_inv*GM1 );
	IF_DIM_GE_3( *d2fds2_ptr[ind++] +=  0.0         );

	// d2fds2[:,3,3,3]
	IF_DIM_GE_3( *d2fds2_ptr[ind++] +=  0.0         );
	IF_DIM_GE_3( *d2fds2_ptr[ind++] +=  0.0         );
	IF_DIM_GE_3( *d2fds2_ptr[ind++] += -rho_inv*GM3 );

	// d2fds2[:,4,3,3]
	IF_DIM_GE_3( *d2fds2_ptr[ind++] += -GM1*rho_inv*u     );
	IF_DIM_GE_3( *d2fds2_ptr[ind++] += -GM1*rho_inv*v     );
	IF_DIM_GE_3( *d2fds2_ptr[ind++] += -GM1*rho_inv*w*3.0 );

	// d2fds2[:,:,4,3]
	// d2fds2[:,0,4,3]
	IF_DIM_GE_3( *d2fds2_ptr[ind++] += 0.0 );
	IF_DIM_GE_3( *d2fds2_ptr[ind++] += 0.0 );
	IF_DIM_GE_3( *d2fds2_ptr[ind++] += 0.0 );

	// d2fds2[:,1,4,3]
	IF_DIM_GE_3( *d2fds2_ptr[ind++] += 0.0 );
	IF_DIM_GE_3( *d2fds2_ptr[ind++] += 0.0 );
	IF_DIM_GE_3( *d2fds2_ptr[ind++] += 0.0 );

	// d2fds2[:,2,4,3]
	IF_DIM_GE_3( *d2fds2_ptr[ind++] += 0.0 );
	IF_DIM_GE_3( *d2fds2_ptr[ind++] += 0.0 );
	IF_DIM_GE_3( *d2fds2_ptr[ind++] += 0.0 );

	// d2fds2[:,3,4,3]
	IF_DIM_GE_3( *d2fds2_ptr[ind++] += 0.0 );
	IF_DIM_GE_3( *d2fds2_ptr[ind++] += 0.0 );
	IF_DIM_GE_3( *d2fds2_ptr[ind++] += 0.0 );

	// d2fds2[:,4,4,3]
	IF_DIM_GE_3( *d2fds2_ptr[ind++] += 0.0           );
	IF_DIM_GE_3( *d2fds2_ptr[ind++] += 0.0           );
	IF_DIM_GE_3( *d2fds2_ptr[ind++] += GAMMA*rho_inv );

	// d2fds2[:,:,:,4]
	// d2fds2[:,:,0,4]
	// d2fds2[:,0,0,4]
	IF_DIM_GE_1( *d2fds2_ptr[ind++] += 0.0 );
	IF_DIM_GE_2( *d2fds2_ptr[ind++] += 0.0 );
	IF_DIM_GE_3( *d2fds2_ptr[ind++] += 0.0 );

	// d2fds2[:,1,0,4]
	IF_DIM_GE_1( *d2fds2_ptr[ind++] += 0.0 );
	IF_DIM_GE_2( *d2fds2_ptr[ind++] += 0.0 );
	IF_DIM_GE_3( *d2fds2_ptr[ind++] += 0.0 );

	// d2fds2[:,2,0,4]
	IF_DIM_GE_2( *d2fds2_ptr[ind++] += 0.0 );
	IF_DIM_GE_2( *d2fds2_ptr[ind++] += 0.0 );
	IF_DIM_GE_3( *d2fds2_ptr[ind++] += 0.0 );

	// d2fds2[:,3,0,4]
	IF_DIM_GE_3( *d2fds2_ptr[ind++] += 0.0 );
	IF_DIM_GE_3( *d2fds2_ptr[ind++] += 0.0 );
	IF_DIM_GE_3( *d2fds2_ptr[ind++] += 0.0 );

	// d2fds2[:,4,0,4]
	IF_DIM_GE_1( *d2fds2_ptr[ind++] += u*dbeta_dW[DIM+1] );
	IF_DIM_GE_2( *d2fds2_ptr[ind++] += v*dbeta_dW[DIM+1] );
	IF_DIM_GE_3( *d2fds2_ptr[ind++] += w*dbeta_dW[DIM+1] );

	// d2fds2[:,:,1,4]
	// d2fds2[:,0,1,4]
	IF_DIM_GE_1( *d2fds2_ptr[ind++] += 0.0 );
	IF_DIM_GE_2( *d2fds2_ptr[ind++] += 0.0 );
	IF_DIM_GE_3( *d2fds2_ptr[ind++] += 0.0 );

	// d2fds2[:,1,1,4]
	IF_DIM_GE_1( *d2fds2_ptr[ind++] += 0.0 );
	IF_DIM_GE_2( *d2fds2_ptr[ind++] += 0.0 );
	IF_DIM_GE_3( *d2fds2_ptr[ind++] += 0.0 );

	// d2fds2[:,2,1,4]
	IF_DIM_GE_2( *d2fds2_ptr[ind++] += 0.0 );
	IF_DIM_GE_2( *d2fds2_ptr[ind++] += 0.0 );
	IF_DIM_GE_3( *d2fds2_ptr[ind++] += 0.0 );

	// d2fds2[:,3,1,4]
	IF_DIM_GE_3( *d2fds2_ptr[ind++] += 0.0 );
	IF_DIM_GE_3( *d2fds2_ptr[ind++] += 0.0 );
	IF_DIM_GE_3( *d2fds2_ptr[ind++] += 0.0 );

	// d2fds2[:,4,1,4]
	IF_DIM_GE_1( *d2fds2_ptr[ind++] += GAMMA*rho_inv );
	IF_DIM_GE_2( *d2fds2_ptr[ind++] += 0.0           );
	IF_DIM_GE_3( *d2fds2_ptr[ind++] += 0.0           );

	// d2fds2[:,:,2,4]
	// d2fds2[:,0,2,4]
	IF_DIM_GE_1( *d2fds2_ptr[ind++] += 0.0 );
	IF_DIM_GE_2( *d2fds2_ptr[ind++] += 0.0 );
	IF_DIM_GE_3( *d2fds2_ptr[ind++] += 0.0 );

	// d2fds2[:,1,2,4]
	IF_DIM_GE_2( *d2fds2_ptr[ind++] += 0.0 );
	IF_DIM_GE_2( *d2fds2_ptr[ind++] += 0.0 );
	IF_DIM_GE_3( *d2fds2_ptr[ind++] += 0.0 );

	// d2fds2[:,2,2,4]
	IF_DIM_GE_2( *d2fds2_ptr[ind++] += 0.0 );
	IF_DIM_GE_2( *d2fds2_ptr[ind++] += 0.0 );
	IF_DIM_GE_3( *d2fds2_ptr[ind++] += 0.0 );

	// d2fds2[:,3,2,4]
	IF_DIM_GE_3( *d2fds2_ptr[ind++] += 0.0 );
	IF_DIM_GE_3( *d2fds2_ptr[ind++] += 0.0 );
	IF_DIM_GE_3( *d2fds2_ptr[ind++] += 0.0 );

	// d2fds2[:,4,2,4]
	IF_DIM_GE_2( *d2fds2_ptr[ind++] += 0.0           );
	IF_DIM_GE_2( *d2fds2_ptr[ind++] += GAMMA*rho_inv );
	IF_DIM_GE_3( *d2fds2_ptr[ind++] += 0.0           );

	// d2fds2[:,:,3,4]
	// d2fds2[:,0,3,4]
	IF_DIM_GE_3( *d2fds2_ptr[ind++] += 0.0 );
	IF_DIM_GE_3( *d2fds2_ptr[ind++] += 0.0 );
	IF_DIM_GE_3( *d2fds2_ptr[ind++] += 0.0 );

	// d2fds2[:,1,3,4]
	IF_DIM_GE_3( *d2fds2_ptr[ind++] += 0.0 );
	IF_DIM_GE_3( *d2fds2_ptr[ind++] += 0.0 );
	IF_DIM_GE_3( *d2fds2_ptr[ind++] += 0.0 );

	// d2fds2[:,2,3,4]
	IF_DIM_GE_3( *d2fds2_ptr[ind++] += 0.0 );
	IF_DIM_GE_3( *d2fds2_ptr[ind++] += 0.0 );
	IF_DIM_GE_3( *d2fds2_ptr[ind++] += 0.0 );

	// d2fds2[:,3,3,4]
	IF_DIM_GE_3( *d2fds2_ptr[ind++] += 0.0 );
	IF_DIM_GE_3( *d2fds2_ptr[ind++] += 0.0 );
	IF_DIM_GE_3( *d2fds2_ptr[ind++] += 0.0 );

	// d2fds2[:,4,3,4]
	IF_DIM_GE_3( *d2fds2_ptr[ind++] += 0.0           );
	IF_DIM_GE_3( *d2fds2_ptr[ind++] += 0.0           );
	IF_DIM_GE_3( *d2fds2_ptr[ind++] += GAMMA*rho_inv );

	// d2fds2[:,:,4,4]
	// d2fds2[:,0,4,4]
	IF_DIM_GE_1( *d2fds2_ptr[ind++] += 0.0 );
	IF_DIM_GE_2( *d2fds2_ptr[ind++] += 0.0 );
	IF_DIM_GE_3( *d2fds2_ptr[ind++] += 0.0 );

	// d2fds2[:,1,4,4]
	IF_DIM_GE_1( *d2fds2_ptr[ind++] += 0.0 );
	IF_DIM_GE_2( *d2fds2_ptr[ind++] += 0.0 );
	IF_DIM_GE_3( *d2fds2_ptr[ind++] += 0.0 );

	// d2fds2[:,2,4,4]
	IF_DIM_GE_2( *d2fds2_ptr[ind++] += 0.0 );
	IF_DIM_GE_2( *d2fds2_ptr[ind++] += 0.0 );
	IF_DIM_GE_3( *d2fds2_ptr[ind++] += 0.0 );

	// d2fds2[:,3,4,4]
	IF_DIM_GE_3( *d2fds2_ptr[ind++] += 0.0 );
	IF_DIM_GE_3( *d2fds2_ptr[ind++] += 0.0 );
	IF_DIM_GE_3( *d2fds2_ptr[ind++] += 0.0 );

	// d2fds2[:,4,4,4]
	IF_DIM_GE_1( *d2fds2_ptr[ind++] += 0.0 );
	IF_DIM_GE_2( *d2fds2_ptr[ind++] += 0.0 );
	IF_DIM_GE_3( *d2fds2_ptr[ind++] += 0.0 );
}
