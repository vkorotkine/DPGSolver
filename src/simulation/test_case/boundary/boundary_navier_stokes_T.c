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
 *  \brief Provides the templated Navier-Stokes boundary condition functions.
 *  \todo clean-up.
 */

#include <assert.h>
#include <math.h>
#include <float.h>

#include "macros.h"
#include "definitions_bc.h"
#include "definitions_core.h"
#include "definitions_physics.h"


#include "def_templates_boundary.h"

#include "def_templates_multiarray.h"
#include "def_templates_vector.h"

#include "def_templates_face_solver.h"

#include "def_templates_math_functions.h"
#include "def_templates_solution.h"

// Static function declarations ************************************************************************************* //

#define NEQ  NEQ_EULER  ///< Number of equations.
#define NVAR NVAR_EULER ///< Number of variables.

struct Exact_Boundary_Data;

/** \brief Pointer to functions computing the exact velocities on the boundary.
 *
 *  \param xyz     xyz-coordinates of the point.
 *  \param eb_data \ref Exact_Boundary_Data.
 */
typedef const Type* (*compute_uvw_ex_fptr)
	(const Real xyz[DIM],
	 const struct Exact_Boundary_Data*const eb_data
	);

/// \brief Container for data relating to the exact values on the domain boundary.
struct Exact_Boundary_Data {
	Real rho,   ///< Exact boundary density.
	     E,     ///< Exact boundary total 'e'nergy.
	     nf_E,  ///< Exact boundary 'n'ormal 'v'iscous 'f'lux for the total 'e'nergy equation.
	     omega; ///< Exact boundary angular velocity.

	compute_uvw_ex_fptr compute_uvw_ex; ///< Pointer to function computing the exact velocities on the boundary.
};

/// \brief Set the required members of the \ref Exact_Boundary_Data container for the input boundary condition type.
static void set_Exact_Boundary_Data
	(struct Exact_Boundary_Data*const eb_data, ///< \ref Exact_Boundary_Data.
	 bool*const need_input,                    ///< Flag for whether the input needs to be read.
	 const int viscous_bc_type                 ///< The type of viscous boundary condition.
	);

/** \brief Version of \ref constructor_Boundary_Value_fptr_T computing members using the exact values for all variables
 *         and interpolated values for all gradients.
 *
 *  The boundary conditions are computed as:
 *  - s_b = 2*s_exact - s_i (such that 0.5*(s_b+s_i) == s_exact);
 *  - g_b = g_i_corrected;
 *  where 's' represents an arbitrary 's'olution variable, and 'g' represents an arbitrary 'g'radient variable. The
 *  corrected boundary gradients are specified such that the normal boundary viscous fluxes are identical when computed
 *  either as:
 *  1. n_dot_F = n (dot) f_viscous(s_l,g_l); or
 *  2. n_dof_F = n (dot) 0.5 * ( f_viscous(s_l,g_l) + f_viscous(s_b,g_b) ).
 *
 *  \warning This boundary condition **is not well-posed** as it has been proven that only NVAR-1 boundary conditions
 *           should be imposed for the Navier-Stokes equations on a no-slip wall (Table 2, Remark 11,
 *           \cite Nordstrom2005). Consequently, it should used only under special circumstances, such as for the
 *           Taylor-Couette case where an additional condition is required due to the domain being periodic.
 */
static void constructor_Boundary_Value_T_navier_stokes_no_slip_all_general
	(struct Boundary_Value_T*const bv,               ///< See brief.
	 const struct Boundary_Value_Input_T*const bv_i, ///< See brief.
	 const struct Exact_Boundary_Data*const eb_data  ///< \ref Exact_Boundary_Data.
	);

/** \brief Version of \ref constructor_Boundary_Value_fptr_T computing members using the exact values for velocity
 *         variables and energy gradient and interpolated values for the remainder.
 *
 *  The boundary solution values are specified such that the exact boundary values are obtained when taking the
 *  average of the internal and boundary (ghost) states. The boundary gradients are specified such that the normal
 *  boundary viscous fluxes are identical when computed either as:
 *  1. n_dot_F = n (dot) f_viscous(s_l,g_l); or
 *  2. n_dof_F = n (dot) 0.5 * ( f_viscous(s_l,g_l) + f_viscous(s_b,g_b) ).
 *
 *  \warning This boundary condition **may not be well-posed** when the normal viscous energy flux is non-zero (i.e. the
 *  non-adiabatic case). See Parsani et. al \cite Parsani2015 for discussion.
 */
static void constructor_Boundary_Value_T_navier_stokes_no_slip_flux_general
	(struct Boundary_Value_T*const bv,               ///< See brief.
	 const struct Boundary_Value_Input_T*const bv_i, ///< See brief.
	 const struct Exact_Boundary_Data*const eb_data  ///< \ref Exact_Boundary_Data.
	);

// Interface functions ********************************************************************************************** //

void constructor_Boundary_Value_T_navier_stokes_no_slip_all_rotating
	(struct Boundary_Value_T* bv, const struct Boundary_Value_Input_T* bv_i, const struct Solver_Face_T* s_face,
	 const struct Simulation* sim)
{
	UNUSED(sim); UNUSED(s_face);

	static bool need_input = true;
	static struct Exact_Boundary_Data eb_data;
	set_Exact_Boundary_Data(&eb_data,&need_input,NO_SLIP_ALL_ROTATING_RHO_E);

	constructor_Boundary_Value_T_navier_stokes_no_slip_all_general(bv,bv_i,&eb_data);
}

void constructor_Boundary_Value_T_navier_stokes_no_slip_flux_adiabatic
	(struct Boundary_Value_T* bv, const struct Boundary_Value_Input_T* bv_i, const struct Solver_Face_T* s_face,
	 const struct Simulation* sim)
{
	UNUSED(sim); UNUSED(s_face);

	static bool need_input = true;
	static struct Exact_Boundary_Data eb_data;
	set_Exact_Boundary_Data(&eb_data,&need_input,DIABATIC_FLUX_CONSTANT_ZERO);

	constructor_Boundary_Value_T_navier_stokes_no_slip_flux_general(bv,bv_i,&eb_data);
}

void constructor_Boundary_Value_T_navier_stokes_no_slip_flux_diabatic
	(struct Boundary_Value_T* bv, const struct Boundary_Value_Input_T* bv_i, const struct Solver_Face_T* s_face,
	 const struct Simulation* sim)
{
	UNUSED(sim); UNUSED(s_face);

	static bool need_input = true;
	static struct Exact_Boundary_Data eb_data;
	set_Exact_Boundary_Data(&eb_data,&need_input,DIABATIC_FLUX_CONSTANT);

	constructor_Boundary_Value_T_navier_stokes_no_slip_flux_general(bv,bv_i,&eb_data);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

/** \brief Read/set the required solution data of the \ref Exact_Boundary_Data container (assumes: zero velocity,
 *         diabatic boundary). */
static void read_and_set_data_diabatic_flux
	(struct Exact_Boundary_Data*const eb_data ///< \ref Exact_Boundary_Data.
	);

/** \brief Read/set the required solution data of the \ref Exact_Boundary_Data container (assumes: angular velocity is
 *         required). */
static void read_and_set_data_no_slip_rotating
	(struct Exact_Boundary_Data*const eb_data ///< \ref Exact_Boundary_Data.
	);

/** \brief Read/set the required solution data of the \ref Exact_Boundary_Data container (assumes: density and
 *         total energy are required). */
static void read_and_set_data_rho_E
	(struct Exact_Boundary_Data*const eb_data ///< \ref Exact_Boundary_Data.
	);

static void set_Exact_Boundary_Data
	(struct Exact_Boundary_Data*const eb_data, bool*const need_input, const int viscous_bc_type)
{
	if (*need_input) {
		*need_input = false;
		switch (viscous_bc_type) {
		case DIABATIC_FLUX_CONSTANT_ZERO: // fallthrough
		case DIABATIC_FLUX_CONSTANT:
			read_and_set_data_diabatic_flux(eb_data);
			break;
		case NO_SLIP_ALL_ROTATING_RHO_E:
			read_and_set_data_no_slip_rotating(eb_data);
			read_and_set_data_rho_E(eb_data);
			break;
		default:
			EXIT_ERROR("Unsupported: %d\n",viscous_bc_type);
			break;
		}
	}
}

static void constructor_Boundary_Value_T_navier_stokes_no_slip_all_general
	(struct Boundary_Value_T*const bv, const struct Boundary_Value_Input_T*const bv_i,
	 const struct Exact_Boundary_Data*const eb_data)
{
	const bool*const c_m = bv_i->compute_member;

	const struct const_Multiarray_T* s_l = bv_i->s;

	const Type*const rho_l  = get_col_const_Multiarray_T(0,s_l),
	          *const E_l    = get_col_const_Multiarray_T(NVAR-1,s_l),
	          *const rhou_l = get_col_const_Multiarray_T(1,s_l);
	IF_DIM_GE_2( const Type*const rhov_l = get_col_const_Multiarray_T(2,s_l) );
	IF_DIM_GE_3( const Type*const rhow_l = get_col_const_Multiarray_T(3,s_l) );

	const ptrdiff_t n_n = s_l->extents[0];

	const struct const_Multiarray_R* xyz = xyz = bv_i->xyz;

	assert(c_m[0] == true);
	struct Multiarray_T* s = constructor_empty_Multiarray_T('C',2,(ptrdiff_t[]){n_n,NVAR}); // keep

	Type*const rho  = get_col_Multiarray_T(0,s),
	    *const E    = get_col_Multiarray_T(NVAR-1,s),
	    *const rhou = get_col_Multiarray_T(1,s);
	IF_DIM_GE_2( Type*const rhov = get_col_Multiarray_T(2,s) );
	IF_DIM_GE_3( Type*const rhow = get_col_Multiarray_T(3,s) );

	const Type rho_ex = eb_data->rho,
	           E_ex   = eb_data->E;

	for (int n = 0; n < n_n; ++n) {
		const Real xyz_n[] = ARRAY_DIM( get_col_const_Multiarray_R(0,xyz)[n],
		                                get_col_const_Multiarray_R(1,xyz)[n],
		                                get_col_const_Multiarray_R(2,xyz)[n] );
		const Type*const uvw_ex = eb_data->compute_uvw_ex(xyz_n,eb_data);

		rho[n]  = -rho_l[n]  + 2.0*rho_ex;
		IF_DIM_GE_1( rhou[n] = -rhou_l[n] + 2.0*rho_ex*uvw_ex[0] );
		IF_DIM_GE_2( rhov[n] = -rhov_l[n] + 2.0*rho_ex*uvw_ex[1] );
		IF_DIM_GE_3( rhow[n] = -rhow_l[n] + 2.0*rho_ex*uvw_ex[2] );
		E[n]    = -E_l[n]    + 2.0*E_ex;
	}
	bv->s = (struct const_Multiarray_T*)s;

	if (c_m[1] == true) {
		struct Multiarray_T*const ds_ds = constructor_zero_Multiarray_T('C',3,(ptrdiff_t[]){n_n,NVAR,NVAR}); // keep

		for (int vr_b = 0; vr_b < NVAR; ++vr_b) {
		for (int vr_i = 0; vr_i < NVAR; ++vr_i) {
			if (!(vr_b == vr_i))
				continue;
			struct Vector_T ds_V = interpret_Multiarray_slice_as_Vector_T(ds_ds,(ptrdiff_t[]){vr_b,vr_i});
			set_to_value_Vector_T(&ds_V,-1.0);
		}}
		bv->ds_ds = (const struct const_Multiarray_T*) ds_ds;
	}

	if (bv_i->g)
		bv->g = constructor_copy_const_Multiarray_T(bv_i->g); // keep

	if (c_m[3] == true && bv->g) {
		struct Multiarray_T*const dg_dg =
			constructor_zero_Multiarray_T('C',5,(ptrdiff_t[]){n_n,NVAR,DIM,NVAR,DIM}); // keep

		for (int vr_b = 0; vr_b < NVAR; ++vr_b) {
		for (int d_b  = 0; d_b  < DIM;  ++d_b)  {
		for (int vr_i = 0; vr_i < NVAR; ++vr_i) {
		for (int d_i  = 0; d_i  < DIM;  ++d_i)  {
			if (!(vr_b == vr_i && d_b == d_i))
				continue;
			struct Vector_T dg_V = interpret_Multiarray_slice_as_Vector_T(dg_dg,(ptrdiff_t[]){vr_b,d_b,vr_i,d_i});
			set_to_value_Vector_T(&dg_V,1.0);
		}}}}
		bv->dg_dg = (const struct const_Multiarray_T*) dg_dg;
	}

	if (c_m[4])
		bv->dg_ds = NULL;

	assert(c_m[5] == false);
}

static void constructor_Boundary_Value_T_navier_stokes_no_slip_flux_general
	(struct Boundary_Value_T*const bv, const struct Boundary_Value_Input_T*const bv_i,
	 const struct Exact_Boundary_Data*const eb_data)
{
	const bool*const c_m = bv_i->compute_member;

	const struct const_Multiarray_T* s_l = bv_i->s;

	const Type*const rho_l         = get_col_const_Multiarray_T(0,s_l),
	          *const E_l           = get_col_const_Multiarray_T(NVAR-1,s_l),
	          *const rhouvw_l[DIM] = ARRAY_DIM( get_col_const_Multiarray_T(1,s_l),
	                                            get_col_const_Multiarray_T(2,s_l),
	                                            get_col_const_Multiarray_T(3,s_l) );

	const ptrdiff_t n_n = s_l->extents[0];

	const struct const_Multiarray_R* xyz = xyz = bv_i->xyz;

	assert(c_m[0] == true);
	{
		struct Multiarray_T* s = constructor_empty_Multiarray_T('C',2,(ptrdiff_t[]){n_n,NVAR}); // keep

		Type*const rho         = get_col_Multiarray_T(0,s),
		    *const E           = get_col_Multiarray_T(NVAR-1,s),
		    *const rhouvw[DIM] = ARRAY_DIM( get_col_Multiarray_T(1,s),
		                                    get_col_Multiarray_T(2,s),
		                                    get_col_Multiarray_T(3,s) );

		for (int n = 0; n < n_n; ++n) {
			const Real xyz_n[] = ARRAY_DIM( get_col_const_Multiarray_R(0,xyz)[n],
			                                get_col_const_Multiarray_R(1,xyz)[n],
			                                get_col_const_Multiarray_R(2,xyz)[n] );
			const Type*const uvw_ex = eb_data->compute_uvw_ex(xyz_n,eb_data);

			rho[n] = rho_l[n];
			E[n]   = E_l[n];
			for (int d = 0; d < DIM; ++d)
				rhouvw[d][n] = -rhouvw_l[d][n] + 2.0*rho[n]*uvw_ex[d];
		}
		bv->s = (struct const_Multiarray_T*)s;
	}

	if (c_m[1] == true) {
		struct Multiarray_T*const ds_ds = constructor_zero_Multiarray_T('C',3,(ptrdiff_t[]){n_n,NVAR,NVAR}); // keep

		for (int vr_b = 0; vr_b < NVAR; ++vr_b) {
		for (int vr_i = 0; vr_i < NVAR; ++vr_i) {
			if (!(vr_b == vr_i))
				continue;
			struct Vector_T ds_V = interpret_Multiarray_slice_as_Vector_T(ds_ds,(ptrdiff_t[]){vr_b,vr_i});
			const Real ds_ds_vr = ( (vr_b == 0 || vr_b == NVAR-1) ? 1.0 : -1.0 );
			set_to_value_Vector_T(&ds_V,ds_ds_vr);
		}}

		Type*const duvw_drho[DIM] = ARRAY_DIM( get_col_Multiarray_T(1+NVAR*0,ds_ds),
		                                       get_col_Multiarray_T(2+NVAR*0,ds_ds),
		                                       get_col_Multiarray_T(3+NVAR*0,ds_ds) );

		for (int n = 0; n < n_n; ++n) {
			const Real xyz_n[] = ARRAY_DIM( get_col_const_Multiarray_R(0,xyz)[n],
			                                get_col_const_Multiarray_R(1,xyz)[n],
			                                get_col_const_Multiarray_R(2,xyz)[n] );
			const Type*const uvw_ex = eb_data->compute_uvw_ex(xyz_n,eb_data);
			for (int d = 0; d < DIM; ++d)
				duvw_drho[d][n] = 2.0*uvw_ex[d];
		}
		bv->ds_ds = (const struct const_Multiarray_T*) ds_ds;
	}

	if (bv_i->g) {
		const struct const_Multiarray_T*const g_l = bv_i->g;
		struct Multiarray_T*const g = constructor_empty_Multiarray_T('C',3,(ptrdiff_t[]){n_n,NVAR,DIM}); // keep

		const Type*const drho_l[DIM] = ARRAY_DIM( get_col_const_Multiarray_T(0+NVAR*0,g_l),
		                                          get_col_const_Multiarray_T(0+NVAR*1,g_l),
		                                          get_col_const_Multiarray_T(0+NVAR*2,g_l) ),
		          *const drhouvw_l[DIM][DIM] = TENSOR_DIM( get_col_const_Multiarray_T(1+NVAR*0,g_l),
		                                                   get_col_const_Multiarray_T(1+NVAR*1,g_l),
		                                                   get_col_const_Multiarray_T(1+NVAR*2,g_l),
		                                                   get_col_const_Multiarray_T(2+NVAR*0,g_l),
		                                                   get_col_const_Multiarray_T(2+NVAR*1,g_l),
		                                                   get_col_const_Multiarray_T(2+NVAR*2,g_l),
		                                                   get_col_const_Multiarray_T(3+NVAR*0,g_l),
		                                                   get_col_const_Multiarray_T(3+NVAR*1,g_l),
		                                                   get_col_const_Multiarray_T(3+NVAR*2,g_l) );

		Type*const drho[DIM] = ARRAY_DIM( get_col_Multiarray_T(0+NVAR*0,g),
		                                  get_col_Multiarray_T(0+NVAR*1,g),
		                                  get_col_Multiarray_T(0+NVAR*2,g) ),
		    *const dE[DIM]   = ARRAY_DIM( get_col_Multiarray_T(NVAR-1+NVAR*0,g),
		                                  get_col_Multiarray_T(NVAR-1+NVAR*1,g),
		                                  get_col_Multiarray_T(NVAR-1+NVAR*2,g) ),
		    *const drhouvw[DIM][DIM] = TENSOR_DIM( get_col_Multiarray_T(1+NVAR*0,g),
		                                           get_col_Multiarray_T(1+NVAR*1,g),
		                                           get_col_Multiarray_T(1+NVAR*2,g),
		                                           get_col_Multiarray_T(2+NVAR*0,g),
		                                           get_col_Multiarray_T(2+NVAR*1,g),
		                                           get_col_Multiarray_T(2+NVAR*2,g),
		                                           get_col_Multiarray_T(3+NVAR*0,g),
		                                           get_col_Multiarray_T(3+NVAR*1,g),
		                                           get_col_Multiarray_T(3+NVAR*2,g) );

		for (int n = 0; n < n_n; ++n) {
			for (int d_b = 0; d_b < DIM; ++d_b) {
				drho[d_b][n] = drho_l[d_b][n];
				dE[d_b][n] = DBL_MIN; // Cause obvious error if used.
				for (int dx = 0; dx < DIM; ++dx)
					drhouvw[d_b][dx][n] = drhouvw_l[d_b][dx][n];
			}
		}
		bv->g = (const struct const_Multiarray_T*) g;
		bv->nf_E = eb_data->nf_E;
		bv->nf_E_provided = true;
	}

	if (c_m[3] == true && bv->g) {
		struct Multiarray_T*const dg_dg =
			constructor_zero_Multiarray_T('C',5,(ptrdiff_t[]){n_n,NVAR,DIM,NVAR,DIM}); // keep

		Type*const ddrho_ddrho[DIM] = ARRAY_DIM( get_col_Multiarray_T(0+NVAR*(0+DIM*(0+NVAR*(0))),dg_dg),
		                                         get_col_Multiarray_T(0+NVAR*(1+DIM*(0+NVAR*(1))),dg_dg),
		                                         get_col_Multiarray_T(0+NVAR*(2+DIM*(0+NVAR*(2))),dg_dg) );

		Type*const ddrhouvw_ddrhouvw[DIM][DIM] =
			TENSOR_DIM( get_col_Multiarray_T(1+NVAR*(0+DIM*(1+NVAR*(0))),dg_dg),
			            get_col_Multiarray_T(1+NVAR*(1+DIM*(1+NVAR*(1))),dg_dg),
			            get_col_Multiarray_T(1+NVAR*(2+DIM*(1+NVAR*(2))),dg_dg),
			            get_col_Multiarray_T(2+NVAR*(0+DIM*(2+NVAR*(0))),dg_dg),
			            get_col_Multiarray_T(2+NVAR*(1+DIM*(2+NVAR*(1))),dg_dg),
			            get_col_Multiarray_T(2+NVAR*(2+DIM*(2+NVAR*(2))),dg_dg),
			            get_col_Multiarray_T(3+NVAR*(0+DIM*(3+NVAR*(0))),dg_dg),
			            get_col_Multiarray_T(3+NVAR*(1+DIM*(3+NVAR*(1))),dg_dg),
			            get_col_Multiarray_T(3+NVAR*(2+DIM*(3+NVAR*(2))),dg_dg) );

		for (int n = 0; n < n_n; ++n) {
			for (int d = 0; d < DIM; ++d) {
				ddrho_ddrho[d][n] = 1.0;
				for (int d2 = 0; d2 < DIM; ++d2)
					ddrhouvw_ddrhouvw[d][d2][n] = 1.0;
//				ddE_ddE[d][n]     = 0.0;
			}
		}
		bv->dg_dg = (const struct const_Multiarray_T*) dg_dg;
	}

	if (c_m[4] == true && bv->g) {
		/** Note that boundary values are only used to compute the numerical flux values. In the case of the
		 *  general no-slip flux boundary condition treated here, the numerical flux for the energy equation is
		 *  constant. Further, as the boundary values for the gradient terms for the other equations only depend
		 *  on internal solution gradients (i.e. not the internal solution), this term is set to NULL and simply
		 *  not used when computing the normal flux Jacobian for the energy equation. */
		bv->dg_ds = NULL;
	}

	assert(c_m[5] == false);
}

// Level 1 ********************************************************************************************************** //

/** \brief Version of \ref compute_uvw_ex_fptr imposing zero velocity.
 *  \return See brief. */
static const Type* compute_uvw_ex_zero
	(const Real xyz[DIM],                           ///< See brief.
	 const struct Exact_Boundary_Data*const eb_data ///< See brief.
	);

/** \brief Version of \ref compute_uvw_ex_fptr imposing velocity for a rotating boundary.
 *  \return See brief. */
static const Type* compute_uvw_ex_rotating
	(const Real xyz[DIM],                           ///< See brief.
	 const struct Exact_Boundary_Data*const eb_data ///< See brief.
	);

static void read_and_set_data_diabatic_flux (struct Exact_Boundary_Data*const eb_data)
{
	const int count_to_find = 1;
	int count_found = 0;

	char line[STRLEN_MAX];

	int diabatic_flux_type = VISCOUS_BC_INVALID;

	FILE* input_file = fopen_input('s',NULL,NULL); // closed
	while (fgets(line,sizeof(line),input_file)) {
		read_skip_convert_const_i(line,"diabatic_flux_type",&diabatic_flux_type,&count_found);
	}
	fclose(input_file);

	if (count_found != count_to_find)
		EXIT_ERROR("Did not find the required number of variables");

	switch (diabatic_flux_type) {
	case DIABATIC_FLUX_CONSTANT_ZERO: // fallthrough
	case DIABATIC_FLUX_CONSTANT:
		eb_data->nf_E = get_normal_flux_Energy();
		break;
	default:
		EXIT_ERROR("Unsupported: %d\n",diabatic_flux_type);
		break;
	}

	eb_data->compute_uvw_ex = compute_uvw_ex_zero;
}

static void read_and_set_data_no_slip_rotating (struct Exact_Boundary_Data*const eb_data)
{
	const int count_to_find = 1;
	int count_found = 0;

	char line[STRLEN_MAX];

	FILE* input_file = fopen_input('s',NULL,NULL); // closed
	while (fgets(line,sizeof(line),input_file)) {
		read_skip_string_count_d("omega",&count_found,line,&eb_data->omega);
	}
	fclose(input_file);

	if (count_found != count_to_find)
		EXIT_ERROR("Did not find the required number of variables");

	eb_data->compute_uvw_ex = compute_uvw_ex_rotating;
}

static void read_and_set_data_rho_E (struct Exact_Boundary_Data*const eb_data)
{
	const int count_to_find = 2;
	int count_found = 0;

	char line[STRLEN_MAX];

	FILE* input_file = fopen_input('s',NULL,NULL); // closed
	while (fgets(line,sizeof(line),input_file)) {
		read_skip_string_count_d("rho_b",&count_found,line,&eb_data->rho);
		read_skip_string_count_d("E_b",  &count_found,line,&eb_data->E);
	}
	fclose(input_file);

	if (count_found != count_to_find)
		EXIT_ERROR("Did not find the required number of variables");
}

// Level 2 ********************************************************************************************************** //

static const Type* compute_uvw_ex_zero (const Real xyz[DIM], const struct Exact_Boundary_Data*const eb_data)
{
	UNUSED(xyz); UNUSED(eb_data);
	static const Type uvw[DIM] = {0.0,};
	return uvw;
}

static const Type* compute_uvw_ex_rotating (const Real xyz[DIM], const struct Exact_Boundary_Data*const eb_data)
{
	assert(DIM >= 2);
	const Real x = xyz[0],
	           y = xyz[1];
	const Real omega = eb_data->omega;

	const Real r  = sqrt(x*x+y*y),
	           th = atan2(y,x),
		     Vt = omega*r;

	static Type uvw[DMAX] = {0.0,};
	uvw[0] = -sin(th)*Vt;
	uvw[1] =  cos(th)*Vt;
	return uvw;
}
