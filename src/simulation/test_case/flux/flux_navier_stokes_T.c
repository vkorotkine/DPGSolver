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
#include <math.h>
#include <stddef.h>

#include "macros.h"
#include "definitions_core.h"
#include "definitions_physics.h"

#include "def_templates_multiarray.h"

#include "def_templates_flux.h"
#include "def_templates_math_functions.h"

// Static function declarations ************************************************************************************* //

/// \todo Make a separate header for DIM-related macro functions.
///\{ \name Macros related to the DIM constant.
#if DIM == 1
	#define SUM_DIM(a,b,c) ((a))
	#define TENSOR_DIM(a0,a1,a2,b0,b1,b2,c0,c1,c2) { {(a0),}, }
#elif DIM == 2
	#define SUM_DIM(a,b,c) ((a)+(b))
	#define TENSOR_DIM(a0,a1,a2,b0,b1,b2,c0,c1,c2) { {(a0),(a1),}, {(b0),(b1),}, }
#elif DIM == 3
	#define SUM_DIM(a,b,c) ((a)+(b)+(c))
	#define TENSOR_DIM(a0,a1,a2,b0,b1,b2,c0,c1,c2) { {(a0),(a1),(a2),}, {(b0),(b1),(b2),}, {(c0),(c1),(c2),}, }
#endif
///\}

#define NEQ  NEQ_EULER  ///< Number of equations.
#define NVAR NVAR_EULER ///< Number of variables.

struct Flux_Data_Navier_Stokes;

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

/** \brief Pointer to functions computing the value of the viscosity.
 *
 *  \param input_path \ref Flux_Input_T::input_path.
 *  \param rho        The density.
 *  \param rhouvw     The xyz momentum components.
 *  \param E          The total energy.
 */
typedef Type (*compute_mu_fptr)
	(const char*const input_path,
	 const Type rho,
	 const Type*const rhouvw,
	 const Type E
	);

/** \brief Pointer to functions computing the Jacobian of the viscosity wrt the solution.
 *
 *  \param input_path \ref Flux_Input_T::input_path.
 *  \param rho        The density.
 *  \param rhouvw     The xyz momentum components.
 *  \param E          The total energy.
 */
typedef const Type* (*compute_dmu_ds_fptr)
	(const char*const input_path,
	 const Type rho,
	 const Type*const rhouvw,
	 const Type E
	);

/** \brief Pointer to functions computing the Jacobian of the viscosity wrt the solution gradients.
 *
 *  \param input_path \ref Flux_Input_T::input_path.
 *  \param rho        The density.
 *  \param rhouvw     The xyz momentum components.
 *  \param E          The total energy.
 */
typedef const Type*const* (*compute_dmu_dg_fptr)
	(const char*const input_path,
	 const Type rho,
	 const Type*const rhouvw,
	 const Type E
	);

/// \brief Container for partial derivatives wrt the solution and gradient terms for a scalar quantity.
struct Partials_Scalar {
	Type d0;               ///< 0th order derivative (i.e. the scalar value).
	const Type* d1s;       ///< 1st order derivative wrt the 's'olution.
	const Type*const* d1g; ///< 1st order derivative wrt the 'g'radients.
};

/// \brief Container for partial derivatives wrt the solution and gradient terms for a vector quantity.
struct Partials_Vector {
	const Type* d0;              ///< 0th order derivative (i.e. the vector values).
	const Type*const* d1s;       ///< 1st order derivative wrt the 's'olution.
	const Type*const*const* d1g; ///< 1st order derivative wrt the 'g'radients.
};

/// \brief Container for partial derivatives wrt the solution and gradient terms for a tensor quantity.
struct Partials_Tensor {
	const Type*const* d0;              ///< 0th order derivative (i.e. the tensor values).
	const Type*const*const* d1s;       ///< 1st order derivative wrt the 's'olution.
	const Type*const*const*const* d1g; ///< 1st order derivative wrt the 'g'radients.
};

/// \brief Container for common data used to compute the fluxes and their Jacobians/Hessians.
struct Flux_Data_Navier_Stokes {
	Type rho, ///< The density.
	     E,   ///< The total energy.
	     Pr;  ///< The Prandtl Number.

/// \todo Name change when done replacing with structs.
	const struct Partials_Scalar* mu_P; ///< Partial derivatives for the viscosity.

//	const struct Partials_Vector* uvw; ///< Partial derivatives for the velocity.


	Type rho_inv, ///< The inverse density.
	     mu;      ///< The viscosity.

	const Type* rhouvw, ///< The momentum components.
	          * uvw,    ///< The velocity components.
	          * drho,   ///< The gradients of the density.
	          * dTs;    ///< The gradients of the 's'caled 'T'emperature (T_s == P/(rho*(\gamma-1))).

	const Type*const* drhouvw, ///< Gradients of the momentum components.
	          *const* duvw,    ///< Gradients of the velocity components.
	          *const* tau;     ///< Viscous stress tensor.


	const char* input_path; ///< \ref Flux_Input_T::input_path.

	compute_dmu_ds_fptr compute_dmu_ds; ///< \ref compute_dmu_ds_fptr.
	compute_dmu_dg_fptr compute_dmu_dg; ///< \ref compute_dmu_dg_fptr.
};

/** \brief Return the pointer to the appropriate \ref compute_Flux_Navier_Stokes_fptr specialization based on the
 *         required members.
 *  \return See brief. */
static compute_Flux_Navier_Stokes_fptr get_compute_Flux_Navier_Stokes_fptr
	(const bool*const c_m ///< \ref Flux_Input_T::compute_member.
	);

/** \brief Return the pointer to the appropriate \ref compute_mu_fptr specialization based on the viscosity type.
 *  \return See brief. */
static compute_mu_fptr get_compute_mu_fptr
	(const char*const input_path ///< \ref Flux_Input_T::input_path.
	);

/** \brief Return the pointer to the appropriate \ref compute_dmu_ds_fptr specialization based on the viscosity type.
 *  \return See brief. */
static compute_dmu_ds_fptr get_compute_dmu_ds_fptr
	(const char*const input_path ///< \ref Flux_Input_T::input_path.
	);

/** \brief Return the pointer to the appropriate \ref compute_dmu_dg_fptr specialization based on the viscosity type.
 *  \return See brief. */
static compute_dmu_dg_fptr get_compute_dmu_dg_fptr
	(const char*const input_path ///< \ref Flux_Input_T::input_path.
	);

/** \brief Return the 'Pr'andtl number for the current test case.
 *  \return See brief. */
static Real compute_Pr
	(const char*const input_path ///< \ref Flux_Input_T::input_path.
	);

/** \brief Return a statically allocated \ref Partials_Scalar for the viscosity.
 *  \return See brief. */
static struct Partials_Scalar compute_mu_P
	(const Type rho,             ///< The density.
	 const Type*const rhouvw,    ///< The momentum components.
	 const Type E,               ///< The total energy.
	 const bool*const c_m,       ///< Array of flags indicating members to be computed.
	 const char*const input_path ///< \ref Flux_Input_T::input_path.
	);

/** \brief Return a statically allocated \ref Partials_Vector for the velocities.
 *  \return See brief. */
static struct Partials_Vector compute_uvw_P
	(const Type rho_inv,      ///< The inverse density.
	 const Type*const rhouvw, ///< The momentum components.
	 const bool*const c_m     ///< Array of flags indicating members to be computed.
	);

/** \brief Return a statically allocated \ref Partials_Tensor for the velocity gradients.
 *  \return See brief. */
static struct Partials_Tensor compute_duvw_P
	(const Type rho_inv,                     ///< The inverse density.
	 const struct Partials_Vector*const uvw, ///< See \ref compute_uvw_P.
	 const Type*const drho,                  ///< The density gradients.
	 const Type*const*const drhouvw,         ///< The momentum gradients.
	 const bool*const c_m                    ///< Array of flags indicating members to be computed.
	);

/** \brief Return a statically allocated variable holding the values of the viscous stress tensor.
 *  \return See brief. */
static const Type*const* compute_tau
	(const Type mu,              ///< The viscosity.
	 const Type*const*const duvw ///< The gradients of the velocity components.
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

	const struct const_Multiarray_T*const g = flux_i->g;
	const Type*const drho_p[]  = ARRAY_DIM( get_col_const_Multiarray_T(0+NVAR*0,g),
	                                        get_col_const_Multiarray_T(0+NVAR*1,g),
	                                        get_col_const_Multiarray_T(0+NVAR*2,g) ),
	          *const dE_p[]    = ARRAY_DIM( get_col_const_Multiarray_T(NVAR-1+NVAR*0,g),
	                                        get_col_const_Multiarray_T(NVAR-1+NVAR*1,g),
	                                        get_col_const_Multiarray_T(NVAR-1+NVAR*2,g) ),
	          *const drhouvw_p[DIM][DIM] = TENSOR_DIM( get_col_const_Multiarray_T(1+NVAR*0,g),
	                                                   get_col_const_Multiarray_T(1+NVAR*1,g),
	                                                   get_col_const_Multiarray_T(1+NVAR*2,g),
	                                                   get_col_const_Multiarray_T(2+NVAR*0,g),
	                                                   get_col_const_Multiarray_T(2+NVAR*1,g),
	                                                   get_col_const_Multiarray_T(2+NVAR*2,g),
	                                                   get_col_const_Multiarray_T(3+NVAR*0,g),
	                                                   get_col_const_Multiarray_T(3+NVAR*1,g),
	                                                   get_col_const_Multiarray_T(3+NVAR*2,g) );

	IF_DIM_GE_1( const Type*const*const drhou_p = drhouvw_p[0] );
	IF_DIM_GE_2( const Type*const*const drhov_p = drhouvw_p[1] );
	IF_DIM_GE_3( const Type*const*const drhow_p = drhouvw_p[2] );

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

	compute_dmu_ds_fptr compute_dmu_ds = NULL;
	Type* dfds_ptr[DIM*NEQ*NVAR] = { NULL };
	if (c_m[1]) {
		compute_dmu_ds = get_compute_dmu_ds_fptr(flux_i->input_path);

		struct Multiarray_T*const dfds = flux->df_ds;
		for (int vr = 0; vr < NVAR; ++vr)  {
		for (int eq = 0; eq < NEQ; ++eq)  {
		for (int d = 0; d < DIM; ++d) {
			const int ind = d+DIM*(eq+NEQ*(vr));
			dfds_ptr[ind] = get_col_Multiarray_T(ind,dfds);
		}}}
	}

	compute_dmu_dg_fptr compute_dmu_dg = NULL;
	Type* dfdg_ptr[DIM*NEQ*NVAR*DIM] = { NULL };
	if (c_m[2]) {
		compute_dmu_dg = get_compute_dmu_dg_fptr(flux_i->input_path);

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

	const Real Pr = compute_Pr(flux_i->input_path);

	const ptrdiff_t n_n = s->extents[0];
	for (ptrdiff_t n = 0; n < n_n; ++n) {
		const Type rho      = rho_p[n],
		           E        = E_p[n],
		           rhouvw[] = ARRAY_DIM( rhouvw_p[0][n], rhouvw_p[1][n], rhouvw_p[2][n] );

		const Type drho[DIM] = ARRAY_DIM( drho_p[0][n],  drho_p[1][n],  drho_p[2][n]  ),
		           dE[DIM]   = ARRAY_DIM( dE_p[0][n],    dE_p[1][n],    dE_p[2][n]    ),
		           drhouvw_2[DIM][DIM] = TENSOR_DIM( drhou_p[0][n], drhou_p[1][n], drhou_p[2][n],
		                                             drhov_p[0][n], drhov_p[1][n], drhov_p[2][n],
		                                             drhow_p[0][n], drhow_p[1][n], drhow_p[2][n] ),
		           *const drhouvw[DIM] = ARRAY_DIM( drhouvw_2[0], drhouvw_2[1], drhouvw_2[2] );

		const Type rho_inv  = 1.0/rho,
		           rho_inv2 = rho_inv*rho_inv;
		const struct Partials_Scalar mu_P   = compute_mu_P(rho,rhouvw,E,c_m,flux_i->input_path);
		const struct Partials_Vector uvw_P  = compute_uvw_P(rho_inv,rhouvw,c_m);
		const struct Partials_Tensor duvw_P = compute_duvw_P(rho_inv,&uvw_P,drho,drhouvw,c_m);
UNUSED(mu_P);
EXIT_UNSUPPORTED;

		const Type *const uvw = uvw_P.d0;

		const Type*const*const duvw = duvw_P.d0;

//		const struct Partials_Scalar mu_p = compute_mu_(rho,rhouvw,E,c_m,flux_i->input_path);


		const Type mu = mu_P.d0;
		const Type*const*const tau = compute_tau(mu,duvw);

		IF_DIM_GE_1( const Type u = uvw[0] );
		IF_DIM_GE_2( const Type v = uvw[1] );
		IF_DIM_GE_3( const Type w = uvw[2] );

		IF_DIM_GE_1( const Type* du = duvw[0] );
		IF_DIM_GE_2( const Type* dv = duvw[1] );
		IF_DIM_GE_3( const Type* dw = duvw[2] );

		const Type dE_o_rho[] = ARRAY_DIM( rho_inv2*(dE[0]*rho-E*drho[0]),
		                                   rho_inv2*(dE[1]*rho-E*drho[1]),
		                                   rho_inv2*(dE[2]*rho-E*drho[2]) );
		const Type dV2[] = ARRAY_DIM( 2.0*SUM_DIM(u*du[0],v*dv[0],w*dw[0]),
		                              2.0*SUM_DIM(u*du[1],v*dv[1],w*dw[1]),
		                              2.0*SUM_DIM(u*du[2],v*dv[2],w*dw[2]) );
		const Type dTs[] = ARRAY_DIM( dE_o_rho[0]-0.5*dV2[0], dE_o_rho[1]-0.5*dV2[1], dE_o_rho[2]-0.5*dV2[2] );

		struct Flux_Data_Navier_Stokes flux_data =
			{ .rho     = rho,
			  .E       = E,
			  .rho_inv = rho_inv,
			  .mu  = mu,
			  .Pr  = Pr,

			  .drho    = drho,
			  .rhouvw  = rhouvw,
			  .uvw     = uvw,
			  .dTs = dTs,

			  .drhouvw = drhouvw,
			  .duvw    = duvw,
			  .tau     = tau,

			  .input_path = flux_i->input_path,
			  .compute_dmu_ds = compute_dmu_ds,
			  .compute_dmu_dg = compute_dmu_dg,
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

/// \brief Set the "viscosity_type" parameter based on the value in the input file.
static void set_viscosity_type
	(const char*const input_path,  ///< \ref Flux_Input_T::input_path.
	 int*const viscosity_type_ptr, ///< Pointer to the variable.
	 bool*const need_input         ///< Pointer to the flag for whether the input is still needed.
	);

/** \brief Version of \ref compute_mu_fptr for constant viscosity.
 *  \return See brief. */
static Type compute_mu_constant
	(const char*const input_path, ///< See brief.
	 const Type rho,              ///< See brief.
	 const Type*const rhouvw,     ///< See brief.
	 const Type E                 ///< See brief.
	);

/** \brief Version of \ref compute_dmu_ds_fptr for constant viscosity.
 *  \return See brief. */
static const Type* compute_dmu_ds_constant
	(const char*const input_path, ///< See brief.
	 const Type rho,              ///< See brief.
	 const Type*const rhouvw,     ///< See brief.
	 const Type E                 ///< See brief.
	);

/** \brief Version of \ref compute_dmu_dg_fptr for constant viscosity.
 *  \return See brief. */
static const Type*const* compute_dmu_dg_constant
	(const char*const input_path, ///< See brief.
	 const Type rho,              ///< See brief.
	 const Type*const rhouvw,     ///< See brief.
	 const Type E                 ///< See brief.
	);

/** \brief Version of \ref compute_mu_fptr using the Sutherland formula (eq. (1.56), \cite Toro2009).
 *  \return See brief. */
static Type compute_mu_sutherland
	(const char*const input_path, ///< See brief.
	 const Type rho,              ///< See brief.
	 const Type*const rhouvw,     ///< See brief.
	 const Type E                 ///< See brief.
	);

/** \brief Return a statically allocated array holding the values of the velocity.
 *  \return See brief. */
static const Type* compute_uvw
	(const Type rho_inv,     ///< The inverse density.
	 const Type*const rhouvw ///< The momentum components.
	);

/** \brief Return a statically allocated array holding the values of the velocity Jacobians wrt the solution.
 *  \return See brief. */
static const Type*const* compute_duvw_ds
	(const Type rho_inv,  ///< \ref Flux_Data_Navier_Stokes::rho_inv.
	 const Type*const uvw ///< \ref Flux_Data_Navier_Stokes::uvw.
	);

/** \brief Return a statically allocated array holding the values of the velocity derivatives wrt xyz.
 *  \return See brief. */
static const Type*const* compute_duvw
	(const Type rho_inv,            ///< The inverse density.
	 const Type*const uvw,          ///< The velocity components.
	 const Type*const drho,         ///< The density gradients.
	 const Type*const*const drhouvw ///< The momentum gradients.
	);

/** \brief Return a statically allocated array holding the values of the velocity gradient Jacobians wrt the solution.
 *  \return See brief. */
static const Type*const*const* compute_dduvw_ds
	(const Type rho_inv,             ///< \ref Flux_Data_Navier_Stokes::rho_inv.
	 const Type*const uvw,           ///< \ref Flux_Data_Navier_Stokes::uvw.
	 const Type*const drho,          ///< \ref Flux_Data_Navier_Stokes::drho.
	 const Type*const*const drhouvw, ///< \ref Flux_Data_Navier_Stokes::drhouvw.
	 const Type*const*const duvw_ds  ///< \ref Flux_Data_Navier_Stokes::duvw_ds.
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

static compute_mu_fptr get_compute_mu_fptr (const char*const input_path)
{
	static int viscosity_type = VISCOSITY_INVALID;
	static bool need_input = true;
	set_viscosity_type(input_path,&viscosity_type,&need_input);

	switch (viscosity_type) {
		case VISCOSITY_CONSTANT:   return compute_mu_constant;                    break;
		case VISCOSITY_SUTHERLAND: return compute_mu_sutherland;                  break;
		default:                   EXIT_ERROR("Unsupported: %d.",viscosity_type); break;
	};
}

static compute_dmu_ds_fptr get_compute_dmu_ds_fptr (const char*const input_path)
{
	static int viscosity_type = VISCOSITY_INVALID;
	static bool need_input = true;
	set_viscosity_type(input_path,&viscosity_type,&need_input);

	switch (viscosity_type) {
		case VISCOSITY_CONSTANT:   return compute_dmu_ds_constant;                break;
		case VISCOSITY_SUTHERLAND: EXIT_ADD_SUPPORT;                              break;
		default:                   EXIT_ERROR("Unsupported: %d.",viscosity_type); break;
	};
}

static compute_dmu_dg_fptr get_compute_dmu_dg_fptr (const char*const input_path)
{
	static int viscosity_type = VISCOSITY_INVALID;
	static bool need_input = true;
	set_viscosity_type(input_path,&viscosity_type,&need_input);

	switch (viscosity_type) {
		case VISCOSITY_CONSTANT:   return compute_dmu_dg_constant;                break;
		case VISCOSITY_SUTHERLAND: EXIT_ADD_SUPPORT;                              break;
		default:                   EXIT_ERROR("Unsupported: %d.",viscosity_type); break;
	};
}

static Real compute_Pr (const char*const input_path)
{
	static Real Pr = 0.0;

	static bool need_input = true;
	if (need_input) {
		need_input = false;

		const int count_to_find = 1;
		int count_found = 0;

		char line[STRLEN_MAX];
		FILE* input_file = fopen_input(input_path,'s',NULL); // closed
		while (fgets(line,sizeof(line),input_file)) {
			read_skip_string_count_d("Pr",&count_found,line,&Pr);
		}
		fclose(input_file);

		if (count_found != count_to_find)
			EXIT_ERROR("Did not find the required number of variables");
	}
	return Pr;
}

static struct Partials_Scalar compute_mu_P
	(const Type rho, const Type*const rhouvw, const Type E, const bool*const c_m, const char*const input_path)
{
	static struct Partials_Scalar ps;

/// \todo Move function/function pointer declarations to appropriate levels.
	compute_mu_fptr compute_mu = get_compute_mu_fptr(input_path);
	ps.d0 = compute_mu(input_path,rho,rhouvw,E);

	if (c_m[1]) {
		compute_dmu_ds_fptr compute_dmu_ds = get_compute_dmu_ds_fptr(input_path);
		ps.d1s = compute_dmu_ds(input_path,rho,rhouvw,E);
	} else {
		ps.d1s = NULL;
	}

	if (c_m[2]) {
		compute_dmu_dg_fptr compute_dmu_dg = get_compute_dmu_dg_fptr(input_path);
		ps.d1g = compute_dmu_dg(input_path,rho,rhouvw,E);
	} else {
		ps.d1g = NULL;
	}

	assert(c_m[3] == false); // ADD_SUPPORT;
	assert(c_m[4] == false); // ADD_SUPPORT;
	assert(c_m[5] == false); // ADD_SUPPORT;

	return ps;
}

static struct Partials_Vector compute_uvw_P
	(const Type rho_inv, const Type*const rhouvw, const bool*const c_m)
{
	static struct Partials_Vector pv;

	pv.d0  = compute_uvw(rho_inv,rhouvw);
	pv.d1s = ( c_m[1] ? compute_duvw_ds(rho_inv,pv.d0) : NULL );
	assert(c_m[2] == false); // ADD_SUPPORT;

	assert(c_m[3] == false); // ADD_SUPPORT;
	assert(c_m[4] == false); // ADD_SUPPORT;
	assert(c_m[5] == false); // ADD_SUPPORT;

	return pv;
}

static struct Partials_Tensor compute_duvw_P
	(const Type rho_inv, const struct Partials_Vector*const uvw, const Type*const drho,
	 const Type*const*const drhouvw, const bool*const c_m)
{
	static struct Partials_Tensor pt;

	pt.d0  = compute_duvw(rho_inv,uvw->d0,drho,drhouvw);
	pt.d1s = ( c_m[1] ? compute_dduvw_ds(rho_inv,uvw->d0,drho,drhouvw,uvw->d1s) : NULL );
	assert(c_m[2] == false); // ADD_SUPPORT;

	assert(c_m[3] == false); // ADD_SUPPORT;
	assert(c_m[4] == false); // ADD_SUPPORT;
	assert(c_m[5] == false); // ADD_SUPPORT;

	return pt;

}

static const Type*const* compute_tau (const Type mu, const Type*const*const duvw)
{
	IF_DIM_GE_1( const Type* du = duvw[0] );
	IF_DIM_GE_2( const Type* dv = duvw[1] );
	IF_DIM_GE_3( const Type* dw = duvw[2] );

	const Type divV = SUM_DIM(du[0],dv[1],dw[2]);

	static Type tau[DIM][DIM] = {{0,}};
	IF_DIM_GE_1( tau[0][0] = mu*2.0*(du[0]-divV/3.0) );
	IF_DIM_GE_2( tau[0][1] = mu*(dv[0]+du[1]) );
	IF_DIM_GE_3( tau[0][2] = mu*(dw[0]+du[2]) );
	IF_DIM_GE_2( tau[1][0] = tau[0][1] );
	IF_DIM_GE_2( tau[1][1] = mu*2.0*(dv[1]-divV/3.0) );
	IF_DIM_GE_3( tau[1][2] = mu*(dw[1]+dv[2]) );
	IF_DIM_GE_3( tau[2][0] = tau[0][2] );
	IF_DIM_GE_3( tau[2][1] = tau[1][2] );
	IF_DIM_GE_3( tau[2][2] = mu*2.0*(dw[2]-divV/3.0) );

	static const Type*const tau_array[DIM] = ARRAY_DIM( tau[0], tau[1], tau[2] );
	return tau_array;
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

static void set_viscosity_type
	(const char*const input_path, int*const viscosity_type_ptr, bool*const need_input)
{
	if (*need_input) {
		*need_input = false;

		const int count_to_find = 1;
		int count_found = 0;

		char line[STRLEN_MAX];
		FILE* input_file = fopen_input(input_path,'s',NULL); // closed
		while (fgets(line,sizeof(line),input_file)) {
			read_skip_convert_i(line,"viscosity_type",viscosity_type_ptr,&count_found);
		}
		fclose(input_file);

		if (count_found != count_to_find)
			EXIT_ERROR("Did not find the required number of variables");
	}
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

static Type compute_mu_constant (const char*const input_path, const Type rho, const Type*const rhouvw, const Type E)
{
	UNUSED(rho); UNUSED(rhouvw); UNUSED(E);
	static Real mu = 0.0;

	static bool need_input = true;
	if (need_input) {
		need_input = false;

		const int count_to_find = 1;
		int count_found = 0;

		char line[STRLEN_MAX];
		FILE* input_file = fopen_input(input_path,'s',NULL); // closed
		while (fgets(line,sizeof(line),input_file)) {
			read_skip_string_count_d("mu",&count_found,line,&mu);
		}
		fclose(input_file);

		if (count_found != count_to_find)
			EXIT_ERROR("Did not find the required number of variables");
	}
	return mu;
}

static const Type* compute_dmu_ds_constant
	(const char*const input_path, const Type rho, const Type*const rhouvw, const Type E)
{
	UNUSED(input_path); UNUSED(rho); UNUSED(rhouvw); UNUSED(E);

	static const Type dmu_ds[NVAR] = { 0, };
	return dmu_ds;
}

static const Type*const* compute_dmu_dg_constant
	(const char*const input_path, const Type rho, const Type*const rhouvw, const Type E)
{
	UNUSED(input_path); UNUSED(rho); UNUSED(rhouvw); UNUSED(E);

	static const Type dmu_dg[NVAR][DIM] = {{0,}};
	static const Type* dmu_dg_array[NVAR] = {NULL,};
	for (int vr = 0; vr < NVAR; ++vr)
		dmu_dg_array[vr] = dmu_dg[vr];
	return dmu_dg_array;
}

static Type compute_mu_sutherland (const char*const input_path, const Type rho, const Type*const rhouvw, const Type E)
{
	static Real r_s = 0.0;

	static bool need_input = true;
	if (need_input) {
		need_input = false;

		const int count_to_find = 1;
		int count_found = 0;

		char line[STRLEN_MAX];
		FILE* input_file = fopen_input(input_path,'s',NULL); // closed
		while (fgets(line,sizeof(line),input_file)) {
			read_skip_string_count_d("r_s",&count_found,line,&r_s);
		}
		fclose(input_file);

		if (count_found != count_to_find)
			EXIT_ERROR("Did not find the required number of variables");
	}

	const Type V2 = compute_V2_from_rhouvw_T(rho,rhouvw),
	           p = GM1*(E-0.5*rho*V2),
		     T = p/(rho*r_s);

	static const Real c1 = 1.46e-6,
	                  c2 = 112;

	return c1/(1+c2/T)*sqrt_T(T);
}

static const Type* compute_uvw (const Type rho_inv, const Type*const rhouvw)
{
	IF_DIM_GE_1( const Type rhou = rhouvw[0] );
	IF_DIM_GE_2( const Type rhov = rhouvw[1] );
	IF_DIM_GE_3( const Type rhow = rhouvw[2] );

	static Type uvw[DIM] = { 0, };
	IF_DIM_GE_1( uvw[0] = rho_inv*rhou; )
	IF_DIM_GE_2( uvw[1] = rho_inv*rhov; )
	IF_DIM_GE_3( uvw[2] = rho_inv*rhow; )

	return uvw;
}

static const Type*const* compute_duvw_ds (const Type rho_inv, const Type*const uvw)
{
	IF_DIM_GE_1( const Type u = uvw[0]; )
	IF_DIM_GE_2( const Type v = uvw[1]; )
	IF_DIM_GE_3( const Type w = uvw[2]; )

	static Type duvw_ds[DIM][NVAR] = {{0,}};

	int d = 0;
	IF_DIM_GE_1( duvw_ds[d][0]      = -rho_inv*u );
	IF_DIM_GE_1( duvw_ds[d][1]      =  rho_inv   );
//	IF_DIM_GE_2( duvw_ds[d][2]      =  0.0       );
//	IF_DIM_GE_3( duvw_ds[d][3]      =  0.0       );
//	IF_DIM_GE_1( duvw_ds[d][NVAR-1] =  0.0       );

	d = 1;
	IF_DIM_GE_2( duvw_ds[d][0]      = -rho_inv*v );
//	IF_DIM_GE_2( duvw_ds[d][1]      =  0.0       );
	IF_DIM_GE_2( duvw_ds[d][2]      =  rho_inv   );
//	IF_DIM_GE_3( duvw_ds[d][3]      =  0.0       );
//	IF_DIM_GE_2( duvw_ds[d][NVAR-1] =  0.0       );

	d = 2;
	IF_DIM_GE_3( duvw_ds[d][0]      = -rho_inv*w );
//	IF_DIM_GE_3( duvw_ds[d][1]      =  0.0       );
//	IF_DIM_GE_3( duvw_ds[d][2]      =  0.0       );
	IF_DIM_GE_3( duvw_ds[d][3]      =  rho_inv   );
//	IF_DIM_GE_3( duvw_ds[d][NVAR-1] =  0.0       );

	static const Type*const duvw_ds1[DIM] = ARRAY_DIM( duvw_ds[0], duvw_ds[1], duvw_ds[2] );
	return duvw_ds1;
}

static const Type*const* compute_duvw
	(const Type rho_inv, const Type*const uvw, const Type*const drho, const Type*const*const drhouvw)
{
	IF_DIM_GE_1( const Type u = uvw[0] );
	IF_DIM_GE_2( const Type v = uvw[1] );
	IF_DIM_GE_3( const Type w = uvw[2] );

	IF_DIM_GE_1( const Type*const drhou = drhouvw[0] );
	IF_DIM_GE_2( const Type*const drhov = drhouvw[1] );
	IF_DIM_GE_3( const Type*const drhow = drhouvw[2] );

	static Type duvw[DIM][DIM] = {{0,}};

	for (int d = 0; d < DIM; ++d) {
		IF_DIM_GE_1( duvw[0][d] = rho_inv*(drhou[d]-drho[d]*u); )
		IF_DIM_GE_2( duvw[1][d] = rho_inv*(drhov[d]-drho[d]*v); )
		IF_DIM_GE_3( duvw[2][d] = rho_inv*(drhow[d]-drho[d]*w); )
	}

	static const Type*const duvw1[DIM] = ARRAY_DIM( duvw[0], duvw[1], duvw[2] );
	return duvw1;
}

static const Type*const*const* compute_dduvw_ds
	(const Type rho_inv, const Type*const uvw, const Type*const drho, const Type*const*const drhouvw,
	 const Type*const*const duvw_ds)
{
	const Type drho_inv_ds[NVAR] = ARRAY_NVAR( -rho_inv*rho_inv, 0.0, 0.0, 0.0, 0.0 );
	IF_DIM_GE_1( const Type u = uvw[0]; )
	IF_DIM_GE_2( const Type v = uvw[1]; )
	IF_DIM_GE_3( const Type w = uvw[2]; )

	IF_DIM_GE_1( const Type*const drhou = drhouvw[0]; )
	IF_DIM_GE_2( const Type*const drhov = drhouvw[1]; )
	IF_DIM_GE_3( const Type*const drhow = drhouvw[2]; )

	IF_DIM_GE_1( const Type*const du_ds = duvw_ds[0]; )
	IF_DIM_GE_2( const Type*const dv_ds = duvw_ds[1]; )
	IF_DIM_GE_3( const Type*const dw_ds = duvw_ds[2]; )

	static Type dduvw_ds[DIM][DIM][NVAR] = {{{0,}}};

	for (int d = 0; d < DIM; ++d) {
	for (int vr = 0; vr < NVAR; ++vr) {
		IF_DIM_GE_1( dduvw_ds[0][d][vr] = drho_inv_ds[vr]*(drhou[d]-drho[d]*u) + rho_inv*(0.0-drho[d]*du_ds[vr]) );
		IF_DIM_GE_2( dduvw_ds[1][d][vr] = drho_inv_ds[vr]*(drhov[d]-drho[d]*v) + rho_inv*(0.0-drho[d]*dv_ds[vr]) );
		IF_DIM_GE_3( dduvw_ds[2][d][vr] = drho_inv_ds[vr]*(drhow[d]-drho[d]*w) + rho_inv*(0.0-drho[d]*dw_ds[vr]) );
	}}

	static const Type*const dduvw_ds2[DIM][DIM] = TENSOR_DIM( dduvw_ds[0][0], dduvw_ds[0][1], dduvw_ds[0][2],
	                                                          dduvw_ds[1][0], dduvw_ds[1][1], dduvw_ds[1][2],
	                                                          dduvw_ds[2][0], dduvw_ds[2][1], dduvw_ds[2][2] );

	static const Type*const*const duvw_ds1[DIM] = ARRAY_DIM( dduvw_ds2[0], dduvw_ds2[1], dduvw_ds2[2] );
	return duvw_ds1;
}

// Level 2 ********************************************************************************************************** //

/** \brief Return a statically allocated array holding the values of the viscous stress tensor Jacobians wrt the
 *         solution.
 *  \return See brief. */
static const Type*const*const* compute_dtau_ds
	(const Type mu,                        ///< \ref Flux_Data_Navier_Stokes::mu.
	 const Type*const*const duvw,          ///< \ref Flux_Data_Navier_Stokes::duvw.
	 const Type*const dmu_ds,              ///< \ref Flux_Data_Navier_Stokes::dmu_ds.
	 const Type*const*const*const dduvw_ds ///< \ref Flux_Data_Navier_Stokes::dduvw_ds.
	);

static void compute_Flux_navier_stokes_0
	(const struct Flux_Data_Navier_Stokes*const flux_data, Type*const f_ptr[DIM*NEQ])
{
	IF_DIM_GE_1( const Type u = flux_data->uvw[0]; )
	IF_DIM_GE_2( const Type v = flux_data->uvw[1]; )
	IF_DIM_GE_3( const Type w = flux_data->uvw[2]; )

	const Type*const dTs = flux_data->dTs;
	const Type*const*const tau = flux_data->tau;

	const Type mu = flux_data->mu,
	           Pr = flux_data->Pr;

	int ind = 0;

	// Note the warning concerning the use of negated fluxes in \ref compute_Flux_T_navier_stokes. Using "-=" below.

	// f[:,0]
	for (int d = 0; d < DIM; ++d)
		*f_ptr[ind++] -= 0.0;

	// f[:,1]
	for (int d = 0; d < DIM; ++d)
		*f_ptr[ind++] -= tau[0][d];

#if DIM >= 2
	// f[:,2]
	for (int d = 0; d < DIM; ++d)
		*f_ptr[ind++] -= tau[1][d];
#endif
#if DIM >= 3
	// f[:,3]
	for (int d = 0; d < DIM; ++d)
		*f_ptr[ind++] -= tau[2][d];
#endif

	// f[:,4]
	for (int d = 0; d < DIM; ++d)
		*f_ptr[ind++] -= GAMMA/Pr*mu*dTs[d] + SUM_DIM( u*tau[d][0],v*tau[d][1],w*tau[d][2] );
}

static void compute_Flux_navier_stokes_1s
	(const struct Flux_Data_Navier_Stokes*const flux_data, Type*const dfds_ptr[DIM*NEQ*NVAR])
{
	const Type rho     = flux_data->rho,
	           E       = flux_data->E,
	           rho_inv = flux_data->rho_inv,
	           mu = flux_data->mu,
	           Pr = flux_data->Pr,

	          *const drho   = flux_data->drho,
	          *const rhouvw = flux_data->rhouvw,
	          *const uvw = flux_data->uvw,
	          *const dTs = flux_data->dTs,

	          *const*const drhouvw = flux_data->drhouvw,
	          *const*const duvw    = flux_data->duvw,
	          *const*const tau     = flux_data->tau;



	IF_DIM_GE_1( const Type u = uvw[0]; )
	IF_DIM_GE_2( const Type v = uvw[1]; )
	IF_DIM_GE_3( const Type w = uvw[2]; )

	const Type*const*const duvw_ds = compute_duvw_ds(rho_inv,uvw);
	const Type*const*const*const dduvw_ds = compute_dduvw_ds(rho_inv,uvw,drho,drhouvw,duvw_ds);

	const Type*const dmu_ds = flux_data->compute_dmu_ds(flux_data->input_path,rho,rhouvw,E);
	const Type*const*const*const dtau_ds = compute_dtau_ds(mu,duvw,dmu_ds,dduvw_ds);
	const Type*const*const dTs_ds = NULL;

	IF_DIM_GE_1( const Type*const du_ds = duvw_ds[0]; )
	IF_DIM_GE_2( const Type*const dv_ds = duvw_ds[1]; )
	IF_DIM_GE_3( const Type*const dw_ds = duvw_ds[2]; )
UNUSED(dduvw_ds);
EXIT_ADD_SUPPORT;

	int ind = 0;

	// Note the warning concerning the use of negated fluxes in \ref compute_Flux_T_navier_stokes. Using "-=" below.

	// dfds[:,:,0:NVAR-1]
	for (int vr = 0; vr < NVAR; ++vr) {
		// dfds[:,0,vr]
		for (int d = 0; d < DIM; ++d)
			*dfds_ptr[ind++] -= 0.0;

		// dfds[:,1,vr]
		for (int d = 0; d < DIM; ++d)
			*dfds_ptr[ind++] -= dtau_ds[0][d][vr];

#if DIM >= 2
		// dfds[:,2,vr]
		for (int d = 0; d < DIM; ++d)
			*dfds_ptr[ind++] -= dtau_ds[1][d][vr];
#endif
#if DIM >= 3
		// dfds[:,3,vr]
		for (int d = 0; d < DIM; ++d)
			*dfds_ptr[ind++] -= dtau_ds[2][d][vr];
#endif

		// dfds[:,4,vr]
		for (int d = 0; d < DIM; ++d) {
			*dfds_ptr[ind++] -= GAMMA/Pr*( dmu_ds[vr]*dTs[d] + mu*dTs_ds[d][vr] )
		                         + SUM_DIM( du_ds[vr]*tau[d][0] + u*dtau_ds[d][0][vr],
		                                    dv_ds[vr]*tau[d][1] + v*dtau_ds[d][1][vr],
		                                    dw_ds[vr]*tau[d][2] + w*dtau_ds[d][2][vr] );
		}
	}
}

static void compute_Flux_navier_stokes_1g
	(const struct Flux_Data_Navier_Stokes*const flux_data, Type*const dfdg_ptr[DIM*NEQ*NVAR*DIM])
{
	EXIT_ADD_SUPPORT; UNUSED(flux_data); UNUSED(dfdg_ptr);

	// Note the warning concerning the use of negated fluxes in \ref compute_Flux_T_navier_stokes. Using "-=" below.
}

// Level 3 ********************************************************************************************************** //

/** \brief Return a statically allocated array holding the values of the Jacobian of the velocity divergence wrt the
 *         solution.
 *  \return See brief. */
static const Type* compute_ddivV_ds
	(const Type*const*const*const dduvw_ds ///< \ref Flux_Data_Navier_Stokes::dduvw_ds.
	);

static const Type*const*const* compute_dtau_ds
	(const Type mu, const Type*const*const duvw, const Type*const dmu_ds, const Type*const*const*const dduvw_ds)
{
	IF_DIM_GE_1( const Type* du = duvw[0] );
	IF_DIM_GE_2( const Type* dv = duvw[1] );
	IF_DIM_GE_3( const Type* dw = duvw[2] );

	IF_DIM_GE_1( const Type*const*const ddu_ds = dduvw_ds[0]; )
	IF_DIM_GE_2( const Type*const*const ddv_ds = dduvw_ds[1]; )
	IF_DIM_GE_3( const Type*const*const ddw_ds = dduvw_ds[2]; )

	const Type divV = SUM_DIM(du[0],dv[1],dw[2]);
	const Type*const ddivV_ds = compute_ddivV_ds(dduvw_ds);

	static Type dtau_ds[DIM][DIM][NVAR] = {{{0,}}};

	for (int vr = 0; vr < NVAR; ++vr) {
		IF_DIM_GE_1( dtau_ds[0][0][vr] = dmu_ds[vr]*2.0*(du[0]-divV/3.0) + mu*2.0*(ddu_ds[0][vr]-ddivV_ds[vr]/3.0) );
		IF_DIM_GE_2( dtau_ds[0][1][vr] = dmu_ds[vr]*(dv[0]+du[1])        + mu*(ddv_ds[0][vr]+ddu_ds[1][vr])        );
		IF_DIM_GE_3( dtau_ds[0][2][vr] = dmu_ds[vr]*(dw[0]+du[2])        + mu*(ddw_ds[0][vr]+ddu_ds[2][vr])        );
		IF_DIM_GE_2( dtau_ds[1][0][vr] = dtau_ds[0][1][vr] );
		IF_DIM_GE_2( dtau_ds[1][1][vr] = dmu_ds[vr]*2.0*(dv[1]-divV/3.0) + mu*2.0*(ddv_ds[1][vr]-ddivV_ds[vr]/3.0) );
		IF_DIM_GE_3( dtau_ds[1][2][vr] = dmu_ds[vr]*(dw[1]+dv[2])        + mu*(ddw_ds[1][vr]+ddv_ds[2][vr]) );
		IF_DIM_GE_3( dtau_ds[2][0][vr] = dtau_ds[0][2][vr] );
		IF_DIM_GE_3( dtau_ds[2][1][vr] = dtau_ds[1][2][vr] );
		IF_DIM_GE_3( dtau_ds[2][2][vr] = dmu_ds[vr]*2.0*(dw[2]-divV/3.0) + mu*2.0*(ddw_ds[2][vr]-ddivV_ds[vr]/3.0) );
	}

	static const Type*const dtau_ds2[DIM][DIM] = TENSOR_DIM( dtau_ds[0][0], dtau_ds[0][1], dtau_ds[0][2],
	                                                         dtau_ds[1][0], dtau_ds[1][1], dtau_ds[1][2],
	                                                         dtau_ds[2][0], dtau_ds[2][1], dtau_ds[2][2] );

	static const Type*const*const dtau_ds1[DIM] = ARRAY_DIM( dtau_ds2[0], dtau_ds2[1], dtau_ds2[2] );
	return dtau_ds1;
}

// Level 4 ********************************************************************************************************** //

static const Type* compute_ddivV_ds (const Type*const*const*const dduvw_ds)
{
	IF_DIM_GE_1( const Type*const*const ddu_ds = dduvw_ds[0] );
	IF_DIM_GE_2( const Type*const*const ddv_ds = dduvw_ds[1] );
	IF_DIM_GE_3( const Type*const*const ddw_ds = dduvw_ds[2] );

	static Type ddivV_ds[NVAR] = {0,};

	for (int vr = 0; vr < NVAR; ++vr)
		ddivV_ds[vr] = SUM_DIM( ddu_ds[0][vr], ddv_ds[1][vr], ddw_ds[2][vr] );

	return ddivV_ds;

}
