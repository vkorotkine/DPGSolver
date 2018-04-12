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

#include "face_solver.h"

#include <assert.h>
#include "gsl/gsl_math.h"

#include "element_solver.h"
#include "volume_solver.h"

#include "multiarray.h"

#include "boundary_advection.h"
#include "boundary_diffusion.h"
#include "boundary_euler.h"
#include "boundary_navier_stokes.h"
#include "const_cast.h"
#include "geometry.h"
#include "simulation.h"
#include "test_case.h"

#include "test_complex_face_solver.h"
#include "complex_multiarray.h"

// Static function declarations ************************************************************************************* //

// Templated functions ********************************************************************************************** //

#include "def_templates_type_d.h"
#include "face_solver_T.c"

// Interface functions ********************************************************************************************** //

void copy_members_r_to_c_Solver_Face
	(struct Solver_Face_c*const s_face, const struct Solver_Face*const s_face_r, const struct Simulation*const sim)
{
	const_cast_ptrdiff(&s_face->ind_dof,s_face_r->ind_dof);
	const_cast_i(&s_face->p_ref,s_face_r->p_ref);
	const_cast_i(&s_face->ml,s_face_r->ml);
	const_cast_c(&s_face->cub_type,s_face_r->cub_type);

	destructor_derived_Solver_Face_c((struct Face*)s_face);

	s_face->nf_coef = constructor_copy_Multiarray_c_Multiarray_d(s_face_r->nf_coef); // destructed

	const_constructor_move_const_Multiarray_d(
		&s_face->xyz_fc,constructor_copy_const_Multiarray_d(s_face_r->xyz_fc)); // destructed
	const_constructor_move_const_Multiarray_d(
		&s_face->xyz_fc_ex_b,constructor_copy_const_Multiarray_d(s_face_r->xyz_fc_ex_b)); // destructed
	const_constructor_move_const_Multiarray_d(
		&s_face->normals_fc,constructor_copy_const_Multiarray_d(s_face_r->normals_fc)); // destructed
	const_constructor_move_const_Multiarray_d(
		&s_face->jacobian_det_fc,constructor_copy_const_Multiarray_d(s_face_r->jacobian_det_fc)); // destructed

	set_function_pointers_face_num_flux_c(s_face,sim);

	s_face->nf_fc = constructor_copy_const_Multiarray_c_Multiarray_d(s_face_r->nf_fc); // destructed
}

bool face_is_conforming (const struct Solver_Face*const s_face)
{
	const struct Face*const face = (struct Face*) s_face;
	if (face->boundary)
		return true;

	const struct Solver_Volume*const s_vol[2] = { (struct Solver_Volume*) face->neigh_info[0].volume,
	                                              (struct Solver_Volume*) face->neigh_info[1].volume, };
	if ((s_vol[0]->p_ref == s_vol[1]->p_ref) && (s_vol[0]->ml == s_vol[1]->ml))
		return true;
	return false;
}

int get_dominant_geom_vol_side_index (const struct Solver_Face*const s_face)
{
	const struct Face*const face = (struct Face*) s_face;
	assert(!face->boundary);
	const struct Solver_Volume*const s_vol[2] = { (struct Solver_Volume*) face->neigh_info[0].volume,
	                                              (struct Solver_Volume*) face->neigh_info[1].volume, };
	if (s_vol[0]->ml != s_vol[1]->ml)
		return ( (s_vol[0]->ml < s_vol[1]->ml) ? 0 : 1 );
	else if (s_vol[0]->p_ref != s_vol[1]->p_ref)
		return ( (s_vol[0]->p_ref < s_vol[1]->p_ref) ? 0 : 1 );
	else
		return 0;
}

int compute_face_geometry_order (const struct Solver_Face*const s_face)
{
	const struct Face*const face = (struct Face*) s_face;
	assert(!face->boundary);
	const struct Solver_Volume*const s_vol[2] = { (struct Solver_Volume*) face->neigh_info[0].volume,
	                                              (struct Solver_Volume*) face->neigh_info[1].volume, };
	return GSL_MIN(s_vol[0]->p_ref,s_vol[1]->p_ref);
}

int compute_face_reference_order (const struct Solver_Face*const s_face)
{
	const struct Face*const face = (struct Face*) s_face;
	const struct Solver_Volume*const s_vol[2] = { (struct Solver_Volume*) face->neigh_info[0].volume,
	                                              (struct Solver_Volume*) face->neigh_info[1].volume, };
	if (face->boundary)
		return s_vol[0]->p_ref;
	else
		return GSL_MAX(s_vol[0]->p_ref,s_vol[1]->p_ref);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

