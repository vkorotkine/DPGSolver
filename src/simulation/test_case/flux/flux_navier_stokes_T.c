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
	Type Pr; ///< The Prandtl Number.

	const struct Partials_Scalar* mu_p;  ///< Partial derivatives for the viscosity.
	const struct Partials_Vector* uvw_p, ///< Partial derivatives for the velocity components.
	                            * dTs_p; ///< Partial derivatives for the scaled temperature gradients.
	const struct Partials_Tensor* tau_p; ///< Partial derivatives for the viscous stress tensor.
};

/** \brief Return the pointer to the appropriate \ref compute_Flux_Navier_Stokes_fptr specialization based on the
 *         required members.
 *  \return See brief. */
static compute_Flux_Navier_Stokes_fptr get_compute_Flux_Navier_Stokes_fptr
	(const bool*const c_m ///< \ref Flux_Input_T::compute_member.
	);

/** \brief Check whether the viscosity is constant for the current test case.
 *  \return `true` if yes; `false` otherwise. */
static bool check_if_mu_is_const
	(const char*const input_path ///< \ref Flux_Input_T::input_path.
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

/** \brief Return the 'Pr'andtl number for the current test case.
 *  \return See brief. */
static Real compute_Pr
	(const char*const input_path ///< \ref Flux_Input_T::input_path.
	);

/** \brief Return a statically allocated \ref Partials_Scalar for the viscosity.
 *  \return See brief. */
static struct Partials_Scalar compute_mu_p
	(const Type rho,             ///< The density.
	 const Type*const rhouvw,    ///< The momentum components.
	 const Type E,               ///< The total energy.
	 const bool*const c_m,       ///< Array of flags indicating members to be computed.
	 const char*const input_path ///< \ref Flux_Input_T::input_path.
	);

/** \brief Return a statically allocated \ref Partials_Vector for the velocities.
 *  \return See brief. */
static struct Partials_Vector compute_uvw_p
	(const Type rho_inv,      ///< The inverse density.
	 const Type*const rhouvw, ///< The momentum components.
	 const bool*const c_m     ///< Array of flags indicating members to be computed.
	);

/** \brief Return a statically allocated \ref Partials_Tensor for the velocity gradients.
 *  \return See brief. */
static struct Partials_Tensor compute_duvw_p
	(const Type rho_inv,                       ///< The inverse density.
	 const struct Partials_Vector*const uvw_p, ///< See \ref compute_uvw_p.
	 const Type*const drho,                    ///< The density gradients.
	 const Type*const*const drhouvw,           ///< The momentum gradients.
	 const bool*const c_m                      ///< Array of flags indicating members to be computed.
	);

/** \brief Return a statically allocated \ref Partials_Tensor for the viscous stress tensor.
 *  \return See brief. */
static struct Partials_Tensor compute_tau_p
	(const struct Partials_Scalar*const mu_p,   ///< See \ref compute_mu_p.
	 const struct Partials_Tensor*const duvw_p, ///< See \ref compute_duvw_p.
	 const bool*const c_m,                      ///< Array of flags indicating members to be computed.
	 const bool mu_is_const                     ///< Flag indicating whether the viscosity is constant.
	);
/** \brief Return a statically allocated \ref Partials_Vector for the scaled Temperature (Ts == P/(rho*(\gamma-1)))
 *         gradients.
 *  \return See brief. */
static struct Partials_Vector compute_dTs_p
	(const Type rho,                            ///< The density.
	 const Type E,                              ///< The total energy.
	 const Type*const drho,                     ///< The density gradients.
	 const Type*const dE,                       ///< The total energy gradients.
	 const struct Partials_Vector*const uvw_p,  ///< See \ref compute_uvw_p.
	 const struct Partials_Tensor*const duvw_p, ///< See \ref compute_duvw_p.
	 const bool*const c_m                       ///< Array of flags indicating members to be computed.
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

	const char*const input_path = flux_i->input_path;
	const bool mu_is_const = check_if_mu_is_const(input_path);

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

	const Real Pr = compute_Pr(input_path);

	const ptrdiff_t n_n = s->extents[0];
	for (ptrdiff_t n = 0; n < n_n; ++n) {
		const Type rho      = rho_p[n],
		           E        = E_p[n],
		           rhouvw[] = ARRAY_DIM( rhouvw_p[0][n], rhouvw_p[1][n], rhouvw_p[2][n] );

		const Type drho[DIM] = ARRAY_DIM( drho_p[0][n], drho_p[1][n], drho_p[2][n] ),
		           dE[DIM]   = ARRAY_DIM( dE_p[0][n],   dE_p[1][n],   dE_p[2][n]   ),
		           drhouvw_2[DIM][DIM] = TENSOR_DIM( drhou_p[0][n], drhou_p[1][n], drhou_p[2][n],
		                                             drhov_p[0][n], drhov_p[1][n], drhov_p[2][n],
		                                             drhow_p[0][n], drhow_p[1][n], drhow_p[2][n] ),
		           *const drhouvw[DIM] = ARRAY_DIM( drhouvw_2[0], drhouvw_2[1], drhouvw_2[2] );

		const Type rho_inv  = 1.0/rho;
		const struct Partials_Scalar mu_p   = compute_mu_p(rho,rhouvw,E,c_m,input_path);
		const struct Partials_Vector uvw_p  = compute_uvw_p(rho_inv,rhouvw,c_m);
		const struct Partials_Tensor duvw_p = compute_duvw_p(rho_inv,&uvw_p,drho,drhouvw,c_m),
		                             tau_p  = compute_tau_p(&mu_p,&duvw_p,c_m,mu_is_const);

		const struct Partials_Vector dTs_p  = compute_dTs_p(rho,E,drho,dE,&uvw_p,&duvw_p,c_m);

		struct Flux_Data_Navier_Stokes flux_data =
			{ .Pr    = Pr,
			  .mu_p  = &mu_p,
			  .uvw_p = &uvw_p,
			  .dTs_p = &dTs_p,
			  .tau_p = &tau_p,
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

/** \brief Return a statically allocated array holding the values of the velocity gradient Jacobians wrt the gradients.
 *  \return See brief. */
static const Type*const*const*const* compute_dduvw_dg
	(const Type rho_inv,  ///< \ref Flux_Data_Navier_Stokes::rho_inv.
	 const Type*const uvw ///< \ref Flux_Data_Navier_Stokes::uvw.
	);

/** \brief Return a statically allocated variable holding the values of the viscous stress tensor.
 *  \return See brief. */
static const Type*const* compute_tau
	(const Type mu,              ///< The viscosity.
	 const Type*const*const duvw ///< The gradients of the velocity components.
	);

/** \brief Return a statically allocated array holding the values of the viscous stress tensor Jacobians wrt the
 *         solution.
 *  \return See brief. */
static const Type*const*const* compute_dtau_ds
	(const struct Partials_Scalar*const mu_p,   ///< See \ref compute_mu_p.
	 const struct Partials_Tensor*const duvw_p, ///< See \ref compute_duvw_p.
	 const bool mu_is_const                     ///< Defined for \ref compute_tau_p.
	);

/** \brief Return a statically allocated array holding the values of the viscous stress tensor Jacobians wrt the
 *         gradient.
 *  \return See brief. */
static const Type*const*const*const* compute_dtau_dg
	(const struct Partials_Scalar*const mu_p,  ///< See \ref compute_mu_p.
	 const struct Partials_Tensor*const duvw_p ///< See \ref compute_duvw_p.
	);

/** \brief Return a statically allocated array holding the values of the scaled temperature gradients.
 *  \return See brief. */
static const Type* compute_dTs
	(const Type rho,             ///< See \ref compute_dTs_p.
	 const Type E,               ///< See \ref compute_dTs_p.
	 const Type*const uvw,       ///< See \ref compute_dTs_p.
	 const Type*const drho,      ///< See \ref compute_dTs_p.
	 const Type*const dE,        ///< See \ref compute_dTs_p.
	 const Type*const*const duvw ///< See \ref compute_dTs_p.
	);

/** \brief Return a statically allocated array holding the values of the scaled temperature gradient Jacobians wrt the
 *         solution.
 *  \return See brief. */
static const Type*const* compute_ddTs_ds
	(const Type rho,                           ///< See \ref compute_dTs_p.
	 const Type E,                             ///< See \ref compute_dTs_p.
	 const Type*const drho,                    ///< See \ref compute_dTs_p.
	 const Type*const dE,                      ///< See \ref compute_dTs_p.
	 const struct Partials_Vector*const uvw_p, ///< See \ref compute_dTs_p.
	 const struct Partials_Tensor*const duvw_p ///< See \ref compute_dTs_p.
	);

/** \brief Return a statically allocated array holding the values of the scaled temperature gradient Jacobians wrt the
 *         gradients.
 *  \return See brief. */
static const Type*const*const* compute_ddTs_dg
	(const Type rho,                           ///< See \ref compute_dTs_p.
	 const Type E,                             ///< See \ref compute_dTs_p.
	 const struct Partials_Vector*const uvw_p, ///< See \ref compute_dTs_p.
	 const struct Partials_Tensor*const duvw_p ///< See \ref compute_dTs_p.
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

static bool check_if_mu_is_const (const char*const input_path)
{
	static int viscosity_type = VISCOSITY_INVALID;
	static bool need_input = true;
	set_viscosity_type(input_path,&viscosity_type,&need_input);

	switch (viscosity_type) {
	case VISCOSITY_CONSTANT:
		return true;
		break;
	case VISCOSITY_SUTHERLAND:
		return false;
		break;
	default:
		EXIT_ERROR("Unsupported: %d.",viscosity_type);
		break;
	};
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

static struct Partials_Scalar compute_mu_p
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

	ps.d1g = NULL; // All zeros.

	assert(c_m[3] == false); // ADD_SUPPORT;
	assert(c_m[4] == false); // ADD_SUPPORT;
	assert(c_m[5] == false); // ADD_SUPPORT;

	return ps;
}

static struct Partials_Vector compute_uvw_p
	(const Type rho_inv, const Type*const rhouvw, const bool*const c_m)
{
	static struct Partials_Vector pv;

	pv.d0  = compute_uvw(rho_inv,rhouvw);
	pv.d1s = ( c_m[1] ? compute_duvw_ds(rho_inv,pv.d0) : NULL );
	pv.d1g = NULL; // All zeroes.

	assert(c_m[3] == false); // ADD_SUPPORT;
	assert(c_m[4] == false); // ADD_SUPPORT;
	assert(c_m[5] == false); // ADD_SUPPORT;

	return pv;
}

static struct Partials_Tensor compute_duvw_p
	(const Type rho_inv, const struct Partials_Vector*const uvw_p, const Type*const drho,
	 const Type*const*const drhouvw, const bool*const c_m)
{
	static struct Partials_Tensor pt;

	pt.d0  = compute_duvw(rho_inv,uvw_p->d0,drho,drhouvw);
	pt.d1s = ( c_m[1] ? compute_dduvw_ds(rho_inv,uvw_p->d0,drho,drhouvw,uvw_p->d1s) : NULL );
	pt.d1g = ( c_m[2] ? compute_dduvw_dg(rho_inv,uvw_p->d0) : NULL );

	assert(c_m[3] == false); // ADD_SUPPORT;
	assert(c_m[4] == false); // ADD_SUPPORT;
	assert(c_m[5] == false); // ADD_SUPPORT;

	return pt;
}

static struct Partials_Tensor compute_tau_p
	(const struct Partials_Scalar*const mu_p, const struct Partials_Tensor*const duvw_p, const bool*const c_m,
	 const bool mu_is_const)
{
	static struct Partials_Tensor pt;

	pt.d0  = compute_tau(mu_p->d0,duvw_p->d0);
	pt.d1s = ( c_m[1] ? compute_dtau_ds(mu_p,duvw_p,mu_is_const) : NULL );
	pt.d1g = ( c_m[2] ? compute_dtau_dg(mu_p,duvw_p) : NULL );

	assert(c_m[3] == false); // ADD_SUPPORT;
	assert(c_m[4] == false); // ADD_SUPPORT;
	assert(c_m[5] == false); // ADD_SUPPORT;

	return pt;
}

static struct Partials_Vector compute_dTs_p
	(const Type rho, const Type E, const Type*const drho, const Type*const dE,
	 const struct Partials_Vector*const uvw_p, const struct Partials_Tensor*const duvw_p, const bool*const c_m)
{
	static struct Partials_Vector pv;

	pv.d0  = compute_dTs(rho,E,uvw_p->d0,drho,dE,duvw_p->d0);
	pv.d1s = ( c_m[1] ? compute_ddTs_ds(rho,E,drho,dE,uvw_p,duvw_p) : NULL );
	pv.d1g = ( c_m[2] ? compute_ddTs_dg(rho,E,uvw_p,duvw_p) : NULL );

	assert(c_m[3] == false); // ADD_SUPPORT;
	assert(c_m[4] == false); // ADD_SUPPORT;
	assert(c_m[5] == false); // ADD_SUPPORT;

	return pv;
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
	static const Type dmu_ds[NVAR] = {0,};
	return dmu_ds;
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
	static Type duvw[DIM][DIM] = {{0,}};

	for (int ds = 0; ds < DIM; ++ds) {
	for (int dx = 0; dx < DIM; ++dx) {
		duvw[ds][dx] = rho_inv*(drhouvw[ds][dx]-drho[dx]*uvw[ds]);
	}}

	static const Type*const duvw1[DIM] = ARRAY_DIM( duvw[0], duvw[1], duvw[2] );
	return duvw1;
}

static const Type*const*const* compute_dduvw_ds
	(const Type rho_inv, const Type*const uvw, const Type*const drho, const Type*const*const drhouvw,
	 const Type*const*const duvw_ds)
{
	const Type drho_inv_ds[NVAR] = ARRAY_VAR( -rho_inv*rho_inv, 0.0, 0.0, 0.0, 0.0 );

	static Type dduvw_ds[DIM][DIM][NVAR] = {{{0,}}};

	for (int ds = 0; ds < DIM; ++ds) {
	for (int dx = 0; dx < DIM; ++dx) {
	for (int vr = 0; vr < NVAR; ++vr) {
		dduvw_ds[ds][dx][vr] = drho_inv_ds[vr]*(drhouvw[ds][dx]-drho[dx]*uvw[ds])
		                     + rho_inv*(0.0-drho[dx]*duvw_ds[ds][vr]);
	}}}

	static const Type*const dduvw_ds2[DIM][DIM] = TENSOR_DIM( dduvw_ds[0][0], dduvw_ds[0][1], dduvw_ds[0][2],
	                                                          dduvw_ds[1][0], dduvw_ds[1][1], dduvw_ds[1][2],
	                                                          dduvw_ds[2][0], dduvw_ds[2][1], dduvw_ds[2][2] );

	static const Type*const*const duvw_ds1[DIM] = ARRAY_DIM( dduvw_ds2[0], dduvw_ds2[1], dduvw_ds2[2] );
	return duvw_ds1;
}

static const Type*const*const*const* compute_dduvw_dg (const Type rho_inv, const Type*const uvw)
{
	static Type dduvw_dg[DIM][DIM][DIM][NVAR];

	for (int ds = 0; ds < DIM; ++ds) {
	for (int dx = 0; dx < DIM; ++dx) {
	for (int dg = 0; dg < DIM; ++dg) {
	for (int vr = 0; vr < NVAR; ++vr) {
		const Type ddrhouvw_dg = ( (vr == ds+1 && dg == dx) ? 1.0 : 0.0 ),
		           ddrho_dg    = ( (vr == 0  && dg == dx) ? 1.0 : 0.0 );
		dduvw_dg[ds][dx][dg][vr] = rho_inv*(ddrhouvw_dg-ddrho_dg*uvw[ds]);
	}}}}

	static const Type* dduvw_dg_3[DIM][DIM][DIM];
	static const Type*const* dduvw_dg_2[DIM][DIM];
	static const Type*const*const* dduvw_dg_1[DIM];

	static bool needs_set = true;
	if (needs_set) {
		needs_set = false;
		for (int ds = 0; ds < DIM; ++ds) {
			for (int dx = 0; dx < DIM; ++dx) {
				for (int dg = 0; dg < DIM; ++dg) {
					dduvw_dg_3[ds][dx][dg] = dduvw_dg[ds][dx][dg];
				}
				dduvw_dg_2[ds][dx] = dduvw_dg_3[ds][dx];
			}
			dduvw_dg_1[ds] = dduvw_dg_2[ds];
		}
	}
	return dduvw_dg_1;
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
	IF_DIM_GE_2( tau[1][1] = mu*2.0*(dv[1]-divV/3.0) );
	IF_DIM_GE_3( tau[1][2] = mu*(dw[1]+dv[2]) );
	IF_DIM_GE_3( tau[2][2] = mu*2.0*(dw[2]-divV/3.0) );
	IF_DIM_GE_2( tau[1][0] = tau[0][1] );
	IF_DIM_GE_3( tau[2][0] = tau[0][2] );
	IF_DIM_GE_3( tau[2][1] = tau[1][2] );

	static const Type*const tau_array[DIM] = ARRAY_DIM( tau[0], tau[1], tau[2] );
	return tau_array;
}

static const Type*const*const* compute_dtau_ds
	(const struct Partials_Scalar*const mu_p, const struct Partials_Tensor*const duvw_p, const bool mu_is_const)
{
	const Type mu           = mu_p->d0;
	const Type*const dmu_ds = mu_p->d1s;

	IF_DIM_GE_1( const Type*const*const ddu_ds = duvw_p->d1s[0]; )
	IF_DIM_GE_2( const Type*const*const ddv_ds = duvw_p->d1s[1]; )
	IF_DIM_GE_3( const Type*const*const ddw_ds = duvw_p->d1s[2]; )

	static Type dtau_ds[DIM][DIM][NVAR] = {{{0,}}};

	for (int vr = 0; vr < NVAR; ++vr) {
		const Type ddivV_ds = SUM_DIM( ddu_ds[0][vr], ddv_ds[1][vr], ddw_ds[2][vr] );
		IF_DIM_GE_1( dtau_ds[0][0][vr] = mu*2.0*(ddu_ds[0][vr]-ddivV_ds/3.0) );
		IF_DIM_GE_2( dtau_ds[0][1][vr] = mu*(ddv_ds[0][vr]+ddu_ds[1][vr])        );
		IF_DIM_GE_3( dtau_ds[0][2][vr] = mu*(ddw_ds[0][vr]+ddu_ds[2][vr])        );
		IF_DIM_GE_2( dtau_ds[1][1][vr] = mu*2.0*(ddv_ds[1][vr]-ddivV_ds/3.0) );
		IF_DIM_GE_3( dtau_ds[1][2][vr] = mu*(ddw_ds[1][vr]+ddv_ds[2][vr])        );
		IF_DIM_GE_3( dtau_ds[2][2][vr] = mu*2.0*(ddw_ds[2][vr]-ddivV_ds/3.0) );
	}

	if (!mu_is_const) {
		IF_DIM_GE_1( const Type* du = duvw_p->d0[0] );
		IF_DIM_GE_2( const Type* dv = duvw_p->d0[1] );
		IF_DIM_GE_3( const Type* dw = duvw_p->d0[2] );
		const Type divV = SUM_DIM(du[0],dv[1],dw[2]);

		for (int vr = 0; vr < NVAR; ++vr) {
			IF_DIM_GE_1( dtau_ds[0][0][vr] += dmu_ds[vr]*2.0*(du[0]-divV/3.0) );
			IF_DIM_GE_2( dtau_ds[0][1][vr] += dmu_ds[vr]*(dv[0]+du[1])        );
			IF_DIM_GE_3( dtau_ds[0][2][vr] += dmu_ds[vr]*(dw[0]+du[2])        );
			IF_DIM_GE_2( dtau_ds[1][1][vr] += dmu_ds[vr]*2.0*(dv[1]-divV/3.0) );
			IF_DIM_GE_3( dtau_ds[1][2][vr] += dmu_ds[vr]*(dw[1]+dv[2])        );
			IF_DIM_GE_3( dtau_ds[2][2][vr] += dmu_ds[vr]*2.0*(dw[2]-divV/3.0) );
		}
	}

	for (int vr = 0; vr < NVAR; ++vr) {
		IF_DIM_GE_2( dtau_ds[1][0][vr] = dtau_ds[0][1][vr] );
		IF_DIM_GE_3( dtau_ds[2][0][vr] = dtau_ds[0][2][vr] );
		IF_DIM_GE_3( dtau_ds[2][1][vr] = dtau_ds[1][2][vr] );
	}

	static const Type*const dtau_ds2[DIM][DIM] = TENSOR_DIM( dtau_ds[0][0], dtau_ds[0][1], dtau_ds[0][2],
	                                                         dtau_ds[1][0], dtau_ds[1][1], dtau_ds[1][2],
	                                                         dtau_ds[2][0], dtau_ds[2][1], dtau_ds[2][2] );

	static const Type*const*const dtau_ds1[DIM] = ARRAY_DIM( dtau_ds2[0], dtau_ds2[1], dtau_ds2[2] );
	return dtau_ds1;
}

static const Type*const*const*const* compute_dtau_dg
	(const struct Partials_Scalar*const mu_p, const struct Partials_Tensor*const duvw_p)
{
	const Type mu = mu_p->d0;
	assert(mu_p->d1g == NULL); // Add support.

	IF_DIM_GE_1( const Type*const*const*const ddu_dg = duvw_p->d1g[0]; )
	IF_DIM_GE_2( const Type*const*const*const ddv_dg = duvw_p->d1g[1]; )
	IF_DIM_GE_3( const Type*const*const*const ddw_dg = duvw_p->d1g[2]; )

	static Type dtau_dg[DIM][DIM][DIM][NVAR];

	for (int dg = 0; dg < DIM; ++dg) {
	for (int vr = 0; vr < NVAR; ++vr) {
		const Type ddivV_dg = SUM_DIM( ddu_dg[0][dg][vr], ddv_dg[1][dg][vr], ddw_dg[2][dg][vr] );
		IF_DIM_GE_1( dtau_dg[0][0][dg][vr] = mu*2.0*(ddu_dg[0][dg][vr]-ddivV_dg/3.0) );
		IF_DIM_GE_2( dtau_dg[0][1][dg][vr] = mu*(ddv_dg[0][dg][vr]+ddu_dg[1][dg][vr])        );
		IF_DIM_GE_3( dtau_dg[0][2][dg][vr] = mu*(ddw_dg[0][dg][vr]+ddu_dg[2][dg][vr])        );
		IF_DIM_GE_2( dtau_dg[1][1][dg][vr] = mu*2.0*(ddv_dg[1][dg][vr]-ddivV_dg/3.0) );
		IF_DIM_GE_3( dtau_dg[1][2][dg][vr] = mu*(ddw_dg[1][dg][vr]+ddv_dg[2][dg][vr])        );
		IF_DIM_GE_3( dtau_dg[2][2][dg][vr] = mu*2.0*(ddw_dg[2][dg][vr]-ddivV_dg/3.0) );
	}}

	for (int dg = 0; dg < DIM; ++dg) {
	for (int vr = 0; vr < NVAR; ++vr) {
		IF_DIM_GE_2( dtau_dg[1][0][dg][vr] = dtau_dg[0][1][dg][vr] );
		IF_DIM_GE_3( dtau_dg[2][0][dg][vr] = dtau_dg[0][2][dg][vr] );
		IF_DIM_GE_3( dtau_dg[2][1][dg][vr] = dtau_dg[1][2][dg][vr] );
	}}


	static const Type* dtau_dg_3[DIM][DIM][DIM];
	static const Type*const* dtau_dg_2[DIM][DIM];
	static const Type*const*const* dtau_dg_1[DIM];

	static bool needs_set = true;
	if (needs_set) {
		needs_set = false;
		for (int ds = 0; ds < DIM; ++ds) {
			for (int dx = 0; dx < DIM; ++dx) {
				for (int dg = 0; dg < DIM; ++dg) {
					dtau_dg_3[ds][dx][dg] = dtau_dg[ds][dx][dg];
				}
				dtau_dg_2[ds][dx] = dtau_dg_3[ds][dx];
			}
			dtau_dg_1[ds] = dtau_dg_2[ds];
		}
	}
	return dtau_dg_1;
}

static const Type* compute_dTs
	(const Type rho, const Type E, const Type*const uvw, const Type*const drho, const Type*const dE,
	 const Type*const*const duvw)
{
	const Type rho_inv  = 1.0/rho;
	const Type rho_inv2 = rho_inv*rho_inv;

	IF_DIM_GE_1( const Type u = uvw[0] );
	IF_DIM_GE_2( const Type v = uvw[1] );
	IF_DIM_GE_3( const Type w = uvw[2] );

	IF_DIM_GE_1( const Type* du = duvw[0] );
	IF_DIM_GE_2( const Type* dv = duvw[1] );
	IF_DIM_GE_3( const Type* dw = duvw[2] );

	static Type dTs[DIM];
	for (int d = 0; d < DIM; ++d) {
		const Type dE_o_rho = rho_inv2*(dE[d]*rho-E*drho[d]),
		           dV2      = 2.0*SUM_DIM( u*du[d],v*dv[d],w*dw[d] );
		dTs[d] = dE_o_rho-0.5*dV2;
	}

	return dTs;
}

static const Type*const* compute_ddTs_ds
	(const Type rho, const Type E, const Type*const drho, const Type*const dE,
	 const struct Partials_Vector*const uvw_p, const struct Partials_Tensor*const duvw_p)
{
	const Type rho_inv  = 1.0/rho;
	const Type rho_inv2 = rho_inv*rho_inv;

	IF_DIM_GE_1( const Type u = uvw_p->d0[0] );
	IF_DIM_GE_2( const Type v = uvw_p->d0[1] );
	IF_DIM_GE_3( const Type w = uvw_p->d0[2] );

	IF_DIM_GE_1( const Type*const du = duvw_p->d0[0] );
	IF_DIM_GE_2( const Type*const dv = duvw_p->d0[1] );
	IF_DIM_GE_3( const Type*const dw = duvw_p->d0[2] );

	IF_DIM_GE_1( const Type*const du_ds = uvw_p->d1s[0] );
	IF_DIM_GE_2( const Type*const dv_ds = uvw_p->d1s[1] );
	IF_DIM_GE_3( const Type*const dw_ds = uvw_p->d1s[2] );

	IF_DIM_GE_1( const Type*const*const ddu_ds = duvw_p->d1s[0] );
	IF_DIM_GE_2( const Type*const*const ddv_ds = duvw_p->d1s[1] );
	IF_DIM_GE_3( const Type*const*const ddw_ds = duvw_p->d1s[2] );

	static Type ddTs_ds[DIM][NVAR];
	for (int d = 0; d < DIM; ++d) {
	for (int vr = 0; vr < NVAR; ++vr) {
		const Type drho_inv2_ds = ( vr == 0      ? -2.0*rho_inv*rho_inv2 : 0.0 ),
		           drho_ds      = ( vr == 0      ? 1.0 : 0.0 ),
		           dE_ds        = ( vr == NVAR-1 ? 1.0 : 0.0 );
		const Type dE_o_rho_ds = drho_inv2_ds*(dE[d]*rho-E*drho[d])
		                       + rho_inv2*(dE[d]*drho_ds-dE_ds*drho[d]),
		           dV2_ds      = 2.0*SUM_DIM( du_ds[vr]*du[d]+u*ddu_ds[d][vr],
		                                      dv_ds[vr]*dv[d]+v*ddv_ds[d][vr],
		                                      dw_ds[vr]*dw[d]+w*ddw_ds[d][vr] );
		ddTs_ds[d][vr] = dE_o_rho_ds-0.5*dV2_ds;
	}}

	static const Type*const ddTs_ds_1[DIM] = ARRAY_DIM( ddTs_ds[0], ddTs_ds[1], ddTs_ds[2] );
	return ddTs_ds_1;
}

static const Type*const*const* compute_ddTs_dg
	(const Type rho, const Type E, const struct Partials_Vector*const uvw_p,
	 const struct Partials_Tensor*const duvw_p)
{
	const Type rho_inv  = 1.0/rho;
	const Type rho_inv2 = rho_inv*rho_inv;

	IF_DIM_GE_1( const Type u = uvw_p->d0[0] );
	IF_DIM_GE_2( const Type v = uvw_p->d0[1] );
	IF_DIM_GE_3( const Type w = uvw_p->d0[2] );

	IF_DIM_GE_1( const Type*const*const*const ddu_dg = duvw_p->d1g[0] );
	IF_DIM_GE_2( const Type*const*const*const ddv_dg = duvw_p->d1g[1] );
	IF_DIM_GE_3( const Type*const*const*const ddw_dg = duvw_p->d1g[2] );

	static Type ddTs_dg[DIM][DIM][NVAR];
	for (int dx = 0; dx < DIM; ++dx) {
	for (int dg = 0; dg < DIM; ++dg) {
	for (int vr = 0; vr < NVAR; ++vr) {
		const Type ddrho_dg = ( (vr == 0      && dx == dg) ? 1.0 : 0.0 ),
		           ddE_dg   = ( (vr == NVAR-1 && dx == dg) ? 1.0 : 0.0 );
		const Type dE_o_rho_dg = rho_inv2*(ddE_dg*rho-E*ddrho_dg),
		           dV2_dg      = 2.0*SUM_DIM( u*ddu_dg[dx][dg][vr], v*ddv_dg[dx][dg][vr], w*ddw_dg[dx][dg][vr] );
		ddTs_dg[dx][dg][vr] = dE_o_rho_dg-0.5*dV2_dg;
	}}}

	static const Type*const ddTs_dg_2[DIM][DIM] = TENSOR_DIM( ddTs_dg[0][0], ddTs_dg[0][1], ddTs_dg[0][2],
	                                                          ddTs_dg[1][0], ddTs_dg[1][1], ddTs_dg[1][2],
	                                                          ddTs_dg[2][0], ddTs_dg[2][1], ddTs_dg[2][2] );
	static const Type*const*const ddTs_dg_1[DIM] = ARRAY_DIM( ddTs_dg_2[0], ddTs_dg_2[1], ddTs_dg_2[2] );
	return ddTs_dg_1;
}

// Level 2 ********************************************************************************************************** //

static void compute_Flux_navier_stokes_0
	(const struct Flux_Data_Navier_Stokes*const flux_data, Type*const f_ptr[DIM*NEQ])
{
	const Type Pr = flux_data->Pr;
	const Type mu = flux_data->mu_p->d0;
	const Type*const uvw = flux_data->uvw_p->d0;
	const Type*const dTs = flux_data->dTs_p->d0;
	const Type*const*const tau = flux_data->tau_p->d0;

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
		*f_ptr[ind++] -= GAMMA/Pr*mu*dTs[d] + SUM_DIM( uvw[0]*tau[d][0],uvw[1]*tau[d][1],uvw[2]*tau[d][2] );
}

static void compute_Flux_navier_stokes_1s
	(const struct Flux_Data_Navier_Stokes*const flux_data, Type*const dfds_ptr[DIM*NEQ*NVAR])
{
	const Type Pr = flux_data->Pr;
	const Type mu           = flux_data->mu_p->d0,
	          *const dmu_ds = flux_data->mu_p->d1s;

	const Type*const uvw           = flux_data->uvw_p->d0,
	          *const*const duvw_ds = flux_data->uvw_p->d1s;

	const Type*const dTs           = flux_data->dTs_p->d0,
	          *const*const ddTs_ds = flux_data->dTs_p->d1s;

	const Type*const*const tau           = flux_data->tau_p->d0,
	          *const*const*const dtau_ds = flux_data->tau_p->d1s;

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
			*dfds_ptr[ind++] -= GAMMA/Pr*( dmu_ds[vr]*dTs[d] + mu*ddTs_ds[d][vr] )
			                  + SUM_DIM( duvw_ds[0][vr]*tau[d][0] + uvw[0]*dtau_ds[d][0][vr],
			                             duvw_ds[1][vr]*tau[d][1] + uvw[1]*dtau_ds[d][1][vr],
			                             duvw_ds[2][vr]*tau[d][2] + uvw[2]*dtau_ds[d][2][vr] );
		}
	}
}

static void compute_Flux_navier_stokes_1g
	(const struct Flux_Data_Navier_Stokes*const flux_data, Type*const dfdg_ptr[DIM*NEQ*NVAR*DIM])
{
	const Type Pr = flux_data->Pr;
	const Type mu                 = flux_data->mu_p->d0,
	          *const*const dmu_dg = flux_data->mu_p->d1g;
	assert(dmu_dg == NULL); // Add appropriate terms if no longer true.

	const Type*const uvw                 = flux_data->uvw_p->d0,
	          *const*const*const duvw_dg = flux_data->uvw_p->d1g;
	assert(duvw_dg == NULL);

	const Type*const*const*const ddTs_dg = flux_data->dTs_p->d1g;

	const Type*const*const*const*const dtau_dg = flux_data->tau_p->d1g;

	int ind = 0;

	// Note the warning concerning the use of negated fluxes in \ref compute_Flux_T_navier_stokes. Using "-=" below.

	// dfdg[:,:,:,dg]
	for (int dg = 0; dg < DIM; ++dg) {
	// dfdg[:,:,0:NVAR-1,dg]
	for (int vr = 0; vr < NVAR; ++vr) {
		// dfds[:,0,vr,dg]
		for (int d = 0; d < DIM; ++d)
			*dfdg_ptr[ind++] -= 0.0;

		// dfds[:,1,vr,dg]
		for (int d = 0; d < DIM; ++d)
			*dfdg_ptr[ind++] -= dtau_dg[0][d][dg][vr];

#if DIM >= 2
		// dfds[:,2,vr,dg]
		for (int d = 0; d < DIM; ++d)
			*dfdg_ptr[ind++] -= dtau_dg[1][d][dg][vr];
#endif
#if DIM >= 3
		// dfds[:,3,vr,dg]
		for (int d = 0; d < DIM; ++d)
			*dfdg_ptr[ind++] -= dtau_dg[2][d][dg][vr];
#endif

		// dfds[:,4,vr,dg]
		for (int d = 0; d < DIM; ++d) {
			// Note dmu_dg, duvw_dg assumed to be 0.0.
			*dfdg_ptr[ind++] -= GAMMA/Pr*mu*ddTs_dg[d][dg][vr]
			                  + SUM_DIM( uvw[0]*dtau_dg[d][0][dg][vr],
			                             uvw[1]*dtau_dg[d][1][dg][vr],
			                             uvw[2]*dtau_dg[d][2][dg][vr] );
		}
	}}
}
