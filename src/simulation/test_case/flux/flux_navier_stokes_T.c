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

/// \brief Container for common data used to compute the fluxes and their Jacobians/Hessians.
struct Flux_Data_Navier_Stokes {
	Type rho, ///< The density.
	     u,   ///< The velocity in the x direction.
	     v,   ///< The velocity in the y direction.
	     w,   ///< The velocity in the z direction.
	     E,   ///< The total energy.
	     mu,  ///< The viscosity.
	     Pr;  ///< The Prandtl Number.

	const Type* dT_s; ///< The gradients of the 's'caled 'T'emperature (T_s == P/(rho*(\gamma-1))).

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

/** \brief Pointer to functions computing the value of the viscosity for the Navier-Stokes flux.
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

/** \brief Return the pointer to the appropriate \ref compute_Flux_Navier_Stokes_fptr specialization based on the
 *         required members.
 *  \return See brief. */
static compute_Flux_Navier_Stokes_fptr get_compute_Flux_Navier_Stokes_fptr
	(const bool*const c_m ///< \ref Flux_Input_T::compute_member.
	);

/** \brief Return the pointer to the appropriate \ref compute_mu_fptr specialization based on the viscosity type.
 *  \return See brief. */
static compute_mu_fptr get_compute_mu_fptr
	(const struct Flux_Input_T*const flux_i ///< \ref Flux_Input_T.
	);

/** \brief Return the 'Pr'andtl number for the current test case.
 *  \return See brief. */
static Real compute_Pr
	(const char*const input_path ///< \ref Flux_Input_T::input_path.
	);

/** \brief Return a statically allocated variable holding the values of the viscous stress tensor.
 *  \return See brief. */
static const Type*const* compute_tau
	(const Type mu,       ///< The viscosity.
	 const Type*const du, ///< The gradient of the u-velocity.
	 const Type*const dv, ///< The gradient of the v-velocity.
	 const Type*const dw  ///< The gradient of the w-velocity.
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
	          *const drhou_p[] = ARRAY_DIM( get_col_const_Multiarray_T(1+NVAR*0,g),
	                                        get_col_const_Multiarray_T(1+NVAR*1,g),
	                                        get_col_const_Multiarray_T(1+NVAR*2,g) ),
	          *const drhov_p[] = ARRAY_DIM( get_col_const_Multiarray_T(2+NVAR*0,g),
	                                        get_col_const_Multiarray_T(2+NVAR*1,g),
	                                        get_col_const_Multiarray_T(2+NVAR*2,g) ),
	          *const drhow_p[] = ARRAY_DIM( get_col_const_Multiarray_T(3+NVAR*0,g),
	                                        get_col_const_Multiarray_T(3+NVAR*1,g),
	                                        get_col_const_Multiarray_T(3+NVAR*2,g) ),
	          *const dE_p[]    = ARRAY_DIM( get_col_const_Multiarray_T(NVAR-1+NVAR*0,g),
	                                        get_col_const_Multiarray_T(NVAR-1+NVAR*1,g),
	                                        get_col_const_Multiarray_T(NVAR-1+NVAR*2,g) );

	const bool* c_m = flux_i->compute_member;

/// \todo Set `compute_flux_navier_stokes_n ` as a member of \ref Flux_Input_T and initialize in the constructor.
	compute_Flux_Navier_Stokes_fptr compute_flux_navier_stokes_n = get_compute_Flux_Navier_Stokes_fptr(c_m);

	compute_mu_fptr compute_mu = get_compute_mu_fptr(flux_i);

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

	const Real Pr = compute_Pr(flux_i->input_path);

	const ptrdiff_t n_n = s->extents[0];
	for (ptrdiff_t n = 0; n < n_n; ++n) {
		const Type rho      = rho_p[n],
			     rhou     = rhouvw_p[0][n],
			     rhov     = rhouvw_p[1][n],
			     rhow     = rhouvw_p[2][n],
		           E        = E_p[n],

		           rho_inv  = 1.0/rho,
		           rho_inv2 = rho_inv*rho_inv,
			     u        = rho_inv*rhou,
			     v        = rho_inv*rhov,
			     w        = rho_inv*rhow;

		const Type drho[]  = ARRAY_DIM( drho_p[0][n],  drho_p[1][n],  drho_p[2][n]  ),
		           drhou[] = ARRAY_DIM( drhou_p[0][n], drhou_p[1][n], drhou_p[2][n] ),
		           drhov[] = ARRAY_DIM( drhov_p[0][n], drhov_p[1][n], drhov_p[2][n] ),
		           drhow[] = ARRAY_DIM( drhow_p[0][n], drhow_p[1][n], drhow_p[2][n] ),
		           dE[]    = ARRAY_DIM( dE_p[0][n],    dE_p[1][n],    dE_p[2][n]    ),

		           du[] = ARRAY_DIM( rho_inv*(drhou[0]-drho[0]*u),
			                       rho_inv*(drhou[1]-drho[1]*u),
						     rho_inv*(drhou[2]-drho[2]*u) ),
		           dv[] = ARRAY_DIM( rho_inv*(drhov[0]-drho[0]*v),
			                       rho_inv*(drhov[1]-drho[1]*v),
						     rho_inv*(drhov[2]-drho[2]*v) ),
		           dw[] = ARRAY_DIM( rho_inv*(drhow[0]-drho[0]*w),
			                       rho_inv*(drhow[1]-drho[1]*w),
						     rho_inv*(drhow[2]-drho[2]*w) );

		const Type rhouvw[] = ARRAY_DIM( rhou, rhov, rhow );
		const Type mu = compute_mu(flux_i->input_path,rho,rhouvw,E);
		const Type*const*const tau = compute_tau(mu,du,dv,dw);

		const Type dE_o_rho[] = ARRAY_DIM( rho_inv2*(dE[0]*rho-E*drho[0]),
		                                   rho_inv2*(dE[1]*rho-E*drho[1]),
		                                   rho_inv2*(dE[2]*rho-E*drho[2]) );
		const Type dV2[] = ARRAY_DIM( 2.0*(u*du[0]+v*dv[0]+w*dw[0]),
		                              2.0*(u*du[1]+v*dv[1]+w*dw[1]),
		                              2.0*(u*du[2]+v*dv[2]+w*dw[2]) );
		const Type dT_s[] = ARRAY_DIM( dE_o_rho[0]-0.5*dV2[0], dE_o_rho[1]-0.5*dV2[1], dE_o_rho[2]-0.5*dV2[2] );

		struct Flux_Data_Navier_Stokes flux_data =
			{ .rho  = rho,
			  .u    = u,
			  .v    = v,
			  .w    = w,
			  .E    = E,
			  .mu   = mu,
			  .Pr   = Pr,
			  .dT_s = dT_s,
			  .tau  = tau,
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

/** \brief Version of \ref compute_mu_fptr for constant viscosity.
 *  \return See brief. */
static Type compute_mu_constant
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

static compute_mu_fptr get_compute_mu_fptr (const struct Flux_Input_T*const flux_i)
{
	static int viscosity_type = VISCOSITY_INVALID;

	static bool need_input = true;
	if (need_input) {
		need_input = false;

		const int count_to_find = 1;
		int count_found = 0;

		char line[STRLEN_MAX];
		FILE* input_file = fopen_input(flux_i->input_path,'s',NULL); // closed
		while (fgets(line,sizeof(line),input_file)) {
			read_skip_convert_i(line,"viscosity_type",&viscosity_type,&count_found);
		}
		fclose(input_file);

		if (count_found != count_to_find)
			EXIT_ERROR("Did not find the required number of variables");
	}

	switch (viscosity_type) {
	case VISCOSITY_CONSTANT:
		return compute_mu_constant;
		break;
	case VISCOSITY_SUTHERLAND:
		return compute_mu_sutherland;
		break;
	default:
		EXIT_ERROR("Unsupported: %d.",viscosity_type);
		break;
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

static const Type*const* compute_tau (const Type mu, const Type*const du, const Type*const dv, const Type*const dw)
{
	MAYBE_UNUSED(dv);
	MAYBE_UNUSED(dw);

	const Type divV = 0.0
	IF_DIM_GE_1(    + du[0] )
	IF_DIM_GE_2(    + dv[1] )
	IF_DIM_GE_3(    + dw[2] );

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

	return (const Type*const*) tau;
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

// Level 2 ********************************************************************************************************** //

static void compute_Flux_navier_stokes_0
	(const struct Flux_Data_Navier_Stokes*const flux_data, Type*const f_ptr[DIM*NEQ])
{
	IF_DIM_GE_1( const Type u = flux_data->u; )
	IF_DIM_GE_2( const Type v = flux_data->v; )
	IF_DIM_GE_3( const Type w = flux_data->w; )

	const Type*const dT_s = flux_data->dT_s;
	const Type*const*const tau = flux_data->tau;

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

	// f[:,4]
	IF_DIM_GE_1( *f_ptr[ind++] += - mu*GAMMA/Pr*dT_s[0] - u*tau[0][0] )
	IF_DIM_GE_2(                                        - v*tau[0][1] )
	IF_DIM_GE_3(                                        - w*tau[0][2] );
	IF_DIM_GE_2( *f_ptr[ind++] += - mu*GAMMA/Pr*dT_s[1] - u*tau[1][0] )
	IF_DIM_GE_2(                                        - v*tau[1][1] )
	IF_DIM_GE_3(                                        - w*tau[1][2] );
	IF_DIM_GE_3( *f_ptr[ind++] += - mu*GAMMA/Pr*dT_s[2] - u*tau[2][0] )
	IF_DIM_GE_3(                                        - v*tau[2][1] )
	IF_DIM_GE_3(                                        - w*tau[2][2] );
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
