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
 *  \brief Provides the templated linear advection boundary condition functions.
 */

#include <assert.h>
#include <math.h>

#include "macros.h"
#include "definitions_mesh.h"


#include "def_templates_boundary.h"
#include "boundary_pde_T.c"

#include "def_templates_multiarray.h"

#include "def_templates_face_solver.h"

#include "def_templates_math_functions.h"

// Static function declarations ************************************************************************************* //

// Interface functions ********************************************************************************************** //

void constructor_Boundary_Value_T_advection_inflow
	(struct Boundary_Value_T* bv, const struct Boundary_Value_Input_T* bv_i, const struct Solver_Face_T* face,
	 const struct Simulation* sim)
{
	UNUSED(face);
	const struct const_Multiarray_d*const xyz = bv_i->xyz;
	const bool* c_m = bv_i->compute_member;

	assert(c_m[0] == true);
	bv->s = constructor_sol_bv(xyz,sim); // keep

	if (c_m[1] == true) {
		const ptrdiff_t n_n  = bv->s->extents[0],
		                n_vr = bv->s->extents[1];
		struct Multiarray_T* ds_ds = constructor_empty_Multiarray_T('C',3,(ptrdiff_t[]){n_n,n_vr,n_vr}); // moved
		set_to_value_Multiarray_T(ds_ds,0.0);
		bv->ds_ds = (const struct const_Multiarray_T*) ds_ds; // keep
	}
	assert(c_m[2] == false);
}

void constructor_Boundary_Value_T_advection_outflow
	(struct Boundary_Value_T* bv, const struct Boundary_Value_Input_T* bv_i, const struct Solver_Face_T* face,
	 const struct Simulation* sim)
{
	UNUSED(face);
	UNUSED(sim);

	const bool* c_m = bv_i->compute_member;

	assert(c_m[0] == true);
	bv->s = constructor_copy_const_Multiarray_T(bv_i->s); // keep

	if (c_m[1] == true) {
		const ptrdiff_t n_n  = bv->s->extents[0],
		                n_vr = bv->s->extents[1];
		struct Multiarray_T* ds_ds = constructor_empty_Multiarray_T('C',3,(ptrdiff_t[]){n_n,n_vr,n_vr}); // moved
		set_to_value_Multiarray_T(ds_ds,1.0);
		bv->ds_ds = (const struct const_Multiarray_T*) ds_ds;
	}
	assert(c_m[2] == false);
}

void constructor_Boundary_Value_T_advection_slipwall
	(struct Boundary_Value_T* bv, const struct Boundary_Value_Input_T* bv_i, const struct Solver_Face_T* face,
	 const struct Simulation* sim)
{
	UNUSED(sim);

	static bool need_input = true;
	static struct Sol_Data__Advection sol_data;
	if (need_input) {
		need_input = false;
		read_data_advection(&sol_data);
	}

	const bool* c_m = bv_i->compute_member;
	assert(c_m[0] == true);

	const struct const_Multiarray_T* sol_l = bv_i->s;
	const Type*const u_l = get_col_const_Multiarray_T(0,sol_l);

	const ptrdiff_t n_n = sol_l->extents[0];
	struct Multiarray_T* sol = constructor_empty_Multiarray_T('C',2,(ptrdiff_t[]){n_n,1}); // moved
	Type*const u = get_col_Multiarray_T(0,sol);

	const Real*const xyz[DIM] = ARRAY_DIM( get_col_const_Multiarray_R(0,bv_i->xyz),
	                                       get_col_const_Multiarray_R(1,bv_i->xyz),
	                                       get_col_const_Multiarray_R(2,bv_i->xyz) );
	const struct const_Multiarray_R* normals = bv_i->normals;
	assert(normals->layout == 'R');

	const double h      = ((struct Face*)face)->h;
	const double exp_bn = 1.0;
UNUSED(exp_bn);
	for (int n = 0; n < n_n; n++) {
		const Real xyz_n[DIM] = ARRAY_DIM(xyz[0][n],xyz[1][n],xyz[2][n]);
		const Real*const data_n = get_row_const_Multiarray_d(n,normals),
		          *const b_adv  = sol_data.compute_b_adv(xyz_n);
		const Real b_l2    = norm_R(DIM,b_adv,"L2"),
		           b_dot_n = dot_R(DIM,b_adv,data_n);

//printf("% .15e % .15e\n",b_dot_n,b_l2);
UNUSED(b_l2); UNUSED(b_dot_n);
//		u[n] = u_l[n]*(1.0-2.0*pow_R(b_dot_n/b_l2,2.0)); // original (use vector: b/norm(b)^2)
//		u[n] = u_l[n]*(1.0-2.0*pow_R(b_dot_n/b_l2,1.0)); // modified (use vector: const*(n+b))
//		u[n] = u_l[n]*(1.0-2.0*(1*sqrt(n)+1)*pow_R(fabs(b_dot_n),exp_bn));
const bool condition = 0;//1||(n == n_n-1);
		u[n] = u_l[n]*(1.0- ( condition ? (1*sqrt(n)+1)*pow_R(h,exp_bn) : 0));
//		u[n] = u_l[n]; // outflow
//		u[n] = u_l[n]-2.0*pow_R(fabs(b_dot_n),exp_bn);
//		u[n] = u_l[n]*(1.0-2.0*pow_R(b_dot_n/b_l2,0.0)); // modified (use vector: n/norm(n)^2 == n)
#if TYPE_RC == TYPE_REAL
//printf("%f % e % e\n",b_l2,b_dot_n,-2.0*pow_R(fabs(b_dot_n),exp_bn));
#endif
	}
//printf("\n");
	bv->s = (struct const_Multiarray_T*)sol; // keep

	if (c_m[1] == true) {
		const ptrdiff_t n_n  = bv->s->extents[0],
		                n_vr = bv->s->extents[1];
		struct Multiarray_T* ds_ds = constructor_empty_Multiarray_T('C',3,(ptrdiff_t[]){n_n,n_vr,n_vr}); // moved

		for (int n = 0; n < n_n; n++) {
			const Real xyz_n[DIM] = ARRAY_DIM(xyz[0][n],xyz[1][n],xyz[2][n]);
			const Real*const data_n = get_row_const_Multiarray_d(n,normals),
				    *const b_adv  = sol_data.compute_b_adv(xyz_n);
			const Real b_l2    = norm_R(DIM,b_adv,"L2"),
				     b_dot_n = dot_R(DIM,b_adv,data_n);

UNUSED(b_l2); UNUSED(b_dot_n);
//			ds_ds->data[n] = 1.0-2.0*pow_R(b_dot_n/b_l2,2.0);
//			ds_ds->data[n] = 1.0-2.0*pow_R(b_dot_n/b_l2,1.0);
//			ds_ds->data[n] = 1.0-2.0*(1*sqrt(n)+1)*pow_R(fabs(b_dot_n),exp_bn);
const bool condition = 0;//1||(n == n_n-1);
			ds_ds->data[n] = 1.0- ( condition ? (1*sqrt(n)+1)*pow_R(h,exp_bn) : 0);
//			ds_ds->data[n] = 1.0;
//			ds_ds->data[n] = 1.0;
//			ds_ds->data[n] = 1.0-2.0*pow_R(b_dot_n/b_l2,0.0);
		}
		bv->ds_ds = (const struct const_Multiarray_T*) ds_ds; // keep
	}
	assert(c_m[2] == false);
}

void constructor_Boundary_Value_T_advection_upwind
	(struct Boundary_Value_T* bv, const struct Boundary_Value_Input_T* bv_i, const struct Solver_Face_T* face,
	 const struct Simulation* sim)
{
	UNUSED(face);

	static bool need_input = true;
	static struct Sol_Data__Advection sol_data;
	if (need_input) {
		need_input = false;
		read_data_advection(&sol_data);
	}

	const bool* c_m = bv_i->compute_member;

	assert(c_m[0] == true);
	struct Multiarray_T*const s = (struct Multiarray_T*) constructor_sol_bv(bv_i->xyz,sim); // moved
	const struct const_Multiarray_T*const s_i = constructor_copy_const_Multiarray_T(bv_i->s); // destructed

	const Real*const xyz[DIM] = ARRAY_DIM( get_col_const_Multiarray_R(0,bv_i->xyz),
	                                       get_col_const_Multiarray_R(1,bv_i->xyz),
	                                       get_col_const_Multiarray_R(2,bv_i->xyz) );
	const struct const_Multiarray_R* normals = bv_i->normals;
	assert(normals->layout == 'R');

	const ptrdiff_t n_n = s->extents[0];
	for (int n = 0; n < n_n; n++) {
		const Real xyz_n[DIM] = ARRAY_DIM(xyz[0][n],xyz[1][n],xyz[2][n]);
		const Real*const data_n = get_row_const_Multiarray_d(n,normals),
			    *const b_adv  = sol_data.compute_b_adv(xyz_n);
		const Real b_dot_n = dot_R(DIM,b_adv,data_n);

		if (b_dot_n <= 0.0)
			; // Inflow (do nothing)
		else
			s->data[n] = s_i->data[n];
	}
	destructor_const_Multiarray_T(s_i);
	bv->s = (struct const_Multiarray_T*)s; // keep

	if (c_m[1] == true) {
		const ptrdiff_t n_n  = bv->s->extents[0],
		                n_vr = bv->s->extents[1];
		struct Multiarray_T* ds_ds = constructor_zero_Multiarray_T('C',3,(ptrdiff_t[]){n_n,n_vr,n_vr}); // moved

		for (int n = 0; n < n_n; n++) {
			const Real xyz_n[DIM] = ARRAY_DIM(xyz[0][n],xyz[1][n],xyz[2][n]);
			const Real*const data_n = get_row_const_Multiarray_d(n,normals),
				    *const b_adv  = sol_data.compute_b_adv(xyz_n);
			const Real b_dot_n = dot_R(DIM,b_adv,data_n);

			if (b_dot_n <= 0.0)
				; // Inflow (do nothing)
			else
				ds_ds->data[n] = 1.0;
		}
		bv->ds_ds = (const struct const_Multiarray_T*) ds_ds; // keep
	}
	assert(c_m[2] == false);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //
