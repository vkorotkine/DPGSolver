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

#include "test_complex_boundary.h"

#include "test_complex_operators.h"
#include "test_complex_solve_dg.h"

#include <assert.h>
#include <stdlib.h>

#include "macros.h"

#include "complex_multiarray.h"

#include "element_solver_dg.h"
#include "face_solver.h"
#include "volume_solver_dg_complex.h"

#include "multiarray_operator.h"
#include "operator.h"
#include "simulation.h"

// Static function declarations ************************************************************************************* //

/** \brief `complex` version of \ref constructor_s_fc_interp.
 *  \return See brief. */
static const struct const_Multiarray_c* constructor_s_fc_interp_c
	(const struct Solver_Face* face, ///< See brief.
	 const struct Simulation* sim,   ///< See brief.
	 const int side_index            ///< See brief.
	);

// Interface functions ********************************************************************************************** //

void constructor_Boundary_Value_Input_c_face_s_fcl_interp
	(struct Boundary_Value_Input_c* bv_i, const struct Solver_Face* face, const struct Simulation* sim)
{
	struct Boundary_Value_Input* bv_i_b = (struct Boundary_Value_Input*) bv_i;

	const int side_index = 0;
	bv_i_b->normals = face->normals_fc;
	bv_i_b->xyz     = face->xyz_fc;
	bv_i->s         = constructor_s_fc_interp_c(face,sim,side_index);
	bv_i->g         = NULL;
}

void destructor_Boundary_Value_Input_c (struct Boundary_Value_Input_c* bv_i)
{
	if (bv_i->s)
		destructor_const_Multiarray_c(bv_i->s);
	if (bv_i->g)
		destructor_const_Multiarray_c(bv_i->g);
}

void constructor_Boundary_Value_c_s_fcl_interp
	(struct Boundary_Value_c* bv, const struct Boundary_Value_Input_c* bv_i, const struct Solver_Face* face,
	 const struct Simulation* sim)
{
	UNUSED(bv_i);
	const int side_index = 1;

	struct Multiarray_c* sol_r_fcr = (struct Multiarray_c*) constructor_s_fc_interp_c(face,sim,side_index);
	permute_Multiarray_c_fc(sol_r_fcr,'R',side_index,face);

	bv->s = (const struct const_Multiarray_c*)sol_r_fcr;
	bv->g = NULL;
}

void destructor_Boundary_Value_c (struct Boundary_Value_c* bv)
{
	if (bv->s)
		destructor_const_Multiarray_c(bv->s);
	if (bv->g)
		destructor_const_Multiarray_c(bv->g);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

static const struct const_Multiarray_c* constructor_s_fc_interp_c
	(const struct Solver_Face* s_face, const struct Simulation* sim, const int side_index)
{
	UNUSED(sim);
	const struct Operator* cv0_vs_fc = get_operator__cv0_vs_fc__rlhs_dg(s_face,side_index);

	struct Solver_Volume* s_vol = (struct Solver_Volume*) ((struct Face*)s_face)->neigh_info[side_index].volume;

	struct Complex_DG_Solver_Volume* c_dg_s_vol = (struct Complex_DG_Solver_Volume*) s_vol;
	const struct const_Multiarray_c* s_coef = (const struct const_Multiarray_c*) c_dg_s_vol->sol_coef;

	return constructor_mm_NN1_Operator_const_Multiarray_c(cv0_vs_fc,s_coef,'C','d',s_coef->order,NULL);
}
