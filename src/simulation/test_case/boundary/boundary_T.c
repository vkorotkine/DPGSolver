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

#include "def_templates_boundary.h"

#include "def_templates_multiarray.h"

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
	const int side_index = 0;
	bv_i->normals = s_face->normals_fc;
	bv_i->xyz     = s_face->xyz_fc;
	bv_i->s       = constructor_s_fc_interp(s_face,sim,side_index);
	bv_i->g       = NULL;
}

void constructor_Boundary_Value_Input_face_sg_fcl_interp_T
	(struct Boundary_Value_Input_T* bv_i, const struct Solver_Face_T* s_face, const struct Simulation* sim)
{
	const int side_index = 0;
	bv_i->normals = s_face->normals_fc;
	bv_i->xyz     = s_face->xyz_fc;
	bv_i->s       = constructor_s_fc_interp(s_face,sim,side_index);
	bv_i->g       = constructor_g_fc_interp(s_face,sim,side_index);
}

void destructor_Boundary_Value_Input_T (struct Boundary_Value_Input_T* bv_i)
{
	if (bv_i->s)
		destructor_const_Multiarray_T(bv_i->s);
	if (bv_i->g)
		destructor_const_Multiarray_T(bv_i->g);
}

void constructor_Boundary_Value_s_fcl_interp_T
	(struct Boundary_Value_T* bv, const struct Boundary_Value_Input_T* bv_i, const struct Solver_Face_T* s_face,
	 const struct Simulation* sim)
{
	UNUSED(bv_i);
	const int side_index = 1;

	struct Multiarray_T* sol_r_fcr = (struct Multiarray_T*) constructor_s_fc_interp(s_face,sim,side_index);
	permute_Multiarray_T_fc(sol_r_fcr,'R',side_index,s_face);

	bv->s     = (const struct const_Multiarray_T*)sol_r_fcr;
	bv->g     = NULL;
	bv->ds_ds = NULL;
}

void destructor_Boundary_Value_T (struct Boundary_Value_T* bv)
{
	if (bv->s)
		destructor_const_Multiarray_T(bv->s);
	if (bv->g)
		destructor_const_Multiarray_T(bv->g);
	if (bv->ds_ds)
		destructor_const_Multiarray_T(bv->ds_ds);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

/** \brief Get the pointer to the appropriate \ref Solver_Element::cv0_vr_fc operator.
 *  \return See brief. */
static const struct Operator* get_operator__cv0_vr_fc
	(const int side_index,              ///< The index of the side of the face under consideration.
	 const struct Solver_Face_T* s_face ///< The current \ref Face.
	);

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
	const struct Operator* cv0_vr_fc = get_operator__cv0_vr_fc(side_index,s_face);

	// sim may be used to store a parameter establishing which type of operator to use for the computation.
	UNUSED(sim);
	const char op_format = 'd';

	struct Solver_Volume_T* s_volume = (struct Solver_Volume_T*) ((struct Face*)s_face)->neigh_info[side_index].volume;

	const struct const_Multiarray_T* g_coef = (const struct const_Multiarray_T*) s_volume->grad_coef;

	return constructor_mm_NN1_Operator_const_Multiarray_T(cv0_vr_fc,g_coef,'C',op_format,g_coef->order,NULL);
}

// Level 1 ********************************************************************************************************** //

static const struct Operator* get_operator__cv0_vr_fc (const int side_index, const struct Solver_Face_T* s_face)
{
	const struct Face* face             = (struct Face*) s_face;
	const struct Volume* vol            = face->neigh_info[side_index].volume;
	const struct Solver_Volume_T* s_vol = (struct Solver_Volume_T*) vol;
	const struct Solver_Element* s_e    = (struct Solver_Element*) vol->element;

	const int ind_lf   = face->neigh_info[side_index].ind_lf,
	          ind_href = face->neigh_info[side_index].ind_href;
	const int p_v = s_vol->p_ref,
	          p_f = s_face->p_ref;

	const int curved = ( (s_face->cub_type == 's') ? 0 : 1 );
	return get_Multiarray_Operator(s_e->cv0_vr_fc[curved],(ptrdiff_t[]){ind_lf,ind_href,0,p_f,p_v});
}
