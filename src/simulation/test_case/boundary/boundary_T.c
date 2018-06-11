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
 */

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>

#include "macros.h"
#include "definitions_bc.h"

#include "def_templates_boundary.h"

#include "def_templates_multiarray.h"
#include "def_templates_vector.h"

#include "def_templates_face_solver.h"
#include "def_templates_volume_solver.h"

#include "def_templates_compute_face_rlhs.h"
#include "def_templates_operators.h"
#include "def_templates_solve_dg.h"

// Static function declarations ************************************************************************************* //

/** \brief Constructor for the solution interpolating from the neighbouring volume to the face cubature nodes.
 *  \return See brief. */
static const struct const_Multiarray_T* constructor_s_fc_interp
	(const struct Solver_Face_T* face, ///< Defined for \ref constructor_Boundary_Value_Input_face_s_fcl_interp_T.
	 const struct Simulation* sim,     ///< Defined for \ref constructor_Boundary_Value_Input_face_s_fcl_interp_T.
	 const int side_index              ///< The index of the side of the face under consideration.
	);

/** \brief Constructor for the solution gradient interpolating from the neighbouring volume to the face cubature nodes.
 *  \return See brief. */
static const struct const_Multiarray_T* constructor_g_fc_interp
	(const struct Solver_Face_T* face, ///< Defined for \ref constructor_Boundary_Value_Input_face_s_fcl_interp_T.
	 const struct Simulation* sim,     ///< Defined for \ref constructor_Boundary_Value_Input_face_s_fcl_interp_T.
	 const int side_index              ///< The index of the side of the face under consideration.
	);

// Interface functions ********************************************************************************************** //

void constructor_Boundary_Value_Input_face_s_fcl_interp_T
	(struct Boundary_Value_Input_T* bv_i, const struct Solver_Face_T* s_face, const struct Simulation* sim)
{
	const struct Face*const face = (struct Face*) s_face;
	bv_i->bc = face->bc;
	bv_i->h  = face->h;
	bv_i->p  = s_face->p_ref;

	const int side_index = 0;
	bv_i->normals     = s_face->normals_fc;
	bv_i->normals_std = NULL;
	if (using_exact_normals_for_boundary() && face->boundary && face->bc > BC_CURVED_START) {
		bv_i->normals     = s_face->normals_fc_exact;
		bv_i->normals_std = s_face->normals_fc;
	}
	bv_i->xyz     = s_face->xyz_fc;
	bv_i->xyz_ex  = s_face->xyz_fc_ex_b;
	bv_i->s       = constructor_s_fc_interp(s_face,sim,side_index); // destructed

	bv_i->jacobian_det_fc = s_face->jacobian_det_fc;
}

void constructor_Boundary_Value_Input_face_sg_fcl_interp_T
	(struct Boundary_Value_Input_T* bv_i, const struct Solver_Face_T* s_face, const struct Simulation* sim)
{
	constructor_Boundary_Value_Input_face_s_fcl_interp_T(bv_i,s_face,sim); // destructed

	const int side_index = 0;
	bv_i->g = constructor_g_fc_interp(s_face,sim,side_index); // destructed
}

void destructor_Boundary_Value_Input_T (struct Boundary_Value_Input_T* bv_i)
{
	destructor_conditional_const_Multiarray_T(bv_i->s);
	destructor_conditional_const_Multiarray_T(bv_i->g);
}

void constructor_Boundary_Value_s_fcl_interp_T
	(struct Boundary_Value_T* bv, const struct Boundary_Value_Input_T* bv_i, const struct Solver_Face_T* s_face,
	 const struct Simulation* sim)
{
	UNUSED(bv_i);
	const int side_index = 1;

	struct Multiarray_T* sol_r_fcr = (struct Multiarray_T*) constructor_s_fc_interp(s_face,sim,side_index); // moved
	permute_Multiarray_T_fc(sol_r_fcr,'R',side_index,s_face);

	bv->s     = (const struct const_Multiarray_T*)sol_r_fcr; // destructed
}

void destructor_Boundary_Value_T (struct Boundary_Value_T* bv)
{
	destructor_conditional_const_Multiarray_T(bv->s);
	destructor_conditional_const_Multiarray_T(bv->g);
	destructor_conditional_const_Multiarray_T(bv->ds_ds);
	destructor_conditional_const_Multiarray_T(bv->dg_dg);
	destructor_conditional_const_Multiarray_T(bv->dg_ds);
}

void constructor_Boundary_Value_T_grad_from_internal
	(struct Boundary_Value_T* bv, const struct Boundary_Value_Input_T* bv_i, const int n_var)
{
	const bool*const c_m = bv_i->compute_member;
	const ptrdiff_t n_n = bv_i->s->extents[0];

	if (bv_i->g)
		bv->g = constructor_copy_const_Multiarray_T(bv_i->g); // keep

	if (c_m[3] == true && bv->g) {
		struct Multiarray_T*const dg_dg =
			constructor_zero_Multiarray_T('C',5,(ptrdiff_t[]){n_n,n_var,DIM,n_var,DIM}); // keep

		for (int vr_b = 0; vr_b < n_var; ++vr_b) {
		for (int d_b  = 0; d_b  < DIM;   ++d_b)  {
		for (int vr_i = 0; vr_i < n_var; ++vr_i) {
		for (int d_i  = 0; d_i  < DIM;   ++d_i)  {
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

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

static const struct const_Multiarray_T* constructor_s_fc_interp
	(const struct Solver_Face_T* s_face, const struct Simulation* sim, const int side_index)
{
	const struct Operator* cv0_vs_fc = get_operator__cv0_vs_fc_T(side_index,s_face);

	// sim may be used to store a parameter establishing which type of operator to use for the computation.
	UNUSED(sim);
	const char op_format = 'd';

	struct Solver_Volume_T* s_volume = (struct Solver_Volume_T*) ((struct Face*)s_face)->neigh_info[side_index].volume;

	const struct const_Multiarray_T* s_coef = (const struct const_Multiarray_T*) s_volume->sol_coef;

	return constructor_mm_NN1_Operator_const_Multiarray_T(cv0_vs_fc,s_coef,'C',op_format,s_coef->order,NULL);
}

static const struct const_Multiarray_T* constructor_g_fc_interp
	(const struct Solver_Face_T* s_face, const struct Simulation* sim, const int side_index)
{
	const struct Operator* cv0_vr_fc = get_operator__cv0_vr_fc_T(side_index,s_face);

	// sim may be used to store a parameter establishing which type of operator to use for the computation.
	UNUSED(sim);
	const char op_format = 'd';

	struct Solver_Volume_T* s_volume = (struct Solver_Volume_T*) ((struct Face*)s_face)->neigh_info[side_index].volume;

	const struct const_Multiarray_T* g_coef = (const struct const_Multiarray_T*) s_volume->grad_coef;

	return constructor_mm_NN1_Operator_const_Multiarray_T(cv0_vr_fc,g_coef,'C',op_format,g_coef->order,NULL);
}

#include "undef_templates_boundary.h"

#include "undef_templates_multiarray.h"
#include "undef_templates_vector.h"

#include "undef_templates_face_solver.h"
#include "undef_templates_volume_solver.h"

#include "undef_templates_compute_face_rlhs.h"
#include "undef_templates_operators.h"
#include "undef_templates_solve_dg.h"
