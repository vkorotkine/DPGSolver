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

#include "geometry_blended.h"

#include "definitions_geometry.h"

#include "matrix.h"
#include "multiarray.h"
#include "vector.h"

#include "element_solver.h"
#include "volume.h"
#include "volume_solver.h"

#include "boundary.h"
#include "math_functions.h"
#include "multiarray_operator.h"
#include "operator.h"
#include "simulation.h"
#include "test_case.h"

// Static function declarations ************************************************************************************* //

/** \brief Compute the radius of the curved boundary as that of the first vertex on the boundary.
 *  \return See brief. */
static double compute_radius_from_xyz_ve
	(const int n_component_radius,            /**< The number of components in xyz_ve to use for the computation of
	                                           *   the radius. */
	 const struct const_Matrix_d*const xyz,   ///< The xyz coordinates in which to search.
	 const struct const_Matrix_d*const xyz_ve ///< The xyz coordinates of the vertices.
	);

// Interface functions ********************************************************************************************** //

#include "def_templates_type_d.h"
#include "geometry_blended_T.c"

constructor_xyz_surface_fptr set_constructor_xyz_surface_fptr (const char*const geom_type, const int geom_prm_type)
{
	constructor_xyz_surface_fptr ptr = NULL;

	if (strcmp(geom_type,"n-cylinder") == 0) {
		if      (geom_prm_type == GEOM_PRM_RADIAL_PROJ) ptr = constructor_xyz_surface_cylinder_radial_proj;
		else if (geom_prm_type == GEOM_PRM_ARC_LENGTH ) ptr = constructor_xyz_surface_cylinder_arc_length;
		else if (geom_prm_type == GEOM_PRM_NORMAL_PROJ) ptr = constructor_xyz_surface_cylinder_normal_proj;
		else                                            EXIT_ERROR("Unsupported: %d\n",geom_prm_type);
	} else {
		EXIT_ERROR("Unsupported: %s\n",geom_type);
	}

	return ptr;
}

const struct const_Matrix_d* constructor_xyz_surface_cylinder_radial_proj
	(const struct Blended_Parametric_Data*const b_p_d)
{
	const struct const_Matrix_d xyz_ve_M = interpret_const_Multiarray_as_Matrix_T(b_p_d->xyz_ve);
	const struct const_Matrix_d*const xyz_ve_b =
		constructor_mm_const_Matrix_d('N','N',1.0,b_p_d->vv0_vv_bX->op_std,&xyz_ve_M,'R'); // destructed

	const double r = compute_radius_from_xyz_ve(2,&xyz_ve_M,xyz_ve_b);

	const ptrdiff_t n_n = xyz_ve_b->ext_0,
	                dim = xyz_ve_b->ext_1;
	assert(dim >= 2);

	struct Matrix_d* xyz_surf = constructor_empty_Matrix_d('R',n_n,dim); // returned
	for (int n = 0; n < n_n; ++n) {
		const double*const xyz_b = get_row_const_Matrix_d(n,xyz_ve_b);
		const double theta = atan2(xyz_b[1],xyz_b[0]);

		double*const xyz_s = get_row_Matrix_d(n,xyz_surf);
		xyz_s[0] = r*cos(theta);
		xyz_s[1] = r*sin(theta);
		if (dim == DMAX)
			xyz_s[2] = xyz_b[2];
	}
	destructor_const_Matrix_d(xyz_ve_b);

	transpose_Matrix_d(xyz_surf,true);
	return (struct const_Matrix_d*) xyz_surf;
}

const struct const_Matrix_d* constructor_xyz_surface_cylinder_arc_length
	(const struct Blended_Parametric_Data*const b_p_d)
{
// test for arc length: Should be the same as \theta parametrization for circle.
//	1. Interpolate xyz_ve to boundary (vv0_vv_bv);
//	2. Compute parametric coordinates of vertices
//	3. Interpolate parametric coordinates to curved boundary nodes (vv0_bv_bgc)
//	4. Compute xyz curved
UNUSED(b_p_d); EXIT_ADD_SUPPORT;
}

const struct const_Matrix_d* constructor_xyz_surface_cylinder_normal_proj
	(const struct Blended_Parametric_Data*const b_p_d)
{
//	1. Interpolate xyz_ve to boundary (vv0_vv_bgc);
//	2. Compute normal vector to surface (pointing outwards)
//	3. Compute xyz curved
UNUSED(b_p_d); EXIT_ADD_SUPPORT;
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

/** \brief Find the index of the first node in the matrix of boundary nodes which is also a vertex node.
 *  \return See brief. */
static ptrdiff_t find_boundary_vertex_index
	(const struct const_Matrix_d*const xyz,   ///< The xyz coordinates in which to search.
	 const struct const_Matrix_d*const xyz_ve ///< The xyz coordinates of the vertices.
	);

static double compute_radius_from_xyz_ve
	(const int n_component_radius, const struct const_Matrix_d*const xyz, const struct const_Matrix_d*const xyz_ve)
{
	const ptrdiff_t ind_b = find_boundary_vertex_index(xyz,xyz_ve);
	return norm_d(n_component_radius,get_row_const_Matrix_d(ind_b,xyz),"L2");
}

// Level 1 ********************************************************************************************************** //

static ptrdiff_t find_boundary_vertex_index
	(const struct const_Matrix_d*const xyz, const struct const_Matrix_d*const xyz_ve)
{
	bool found = false;
	ptrdiff_t ind_ve = -1;

	const ptrdiff_t n_n  = xyz->ext_0,
	                n_ve = xyz_ve->ext_0,
	                dim  = xyz->ext_1;
	assert(dim == xyz_ve->ext_1);

	const bool transpose = (xyz_ve->layout != 'R');
	if (transpose)
		transpose_Matrix_d((struct Matrix_d*)xyz_ve,true);

	for (int n = 0; !found && n < n_n; ++n) {
		const double*const data_xyz = get_row_const_Matrix_d(n,xyz);
		for (int ve = 0; ve < n_ve; ++ve) {
			const double*const data_xyz_ve = get_row_const_Matrix_d(ve,xyz_ve);
			const double diff = norm_diff_d(dim,data_xyz_ve,data_xyz,"Inf");
			if (diff < EPS) {
				found = true;
				ind_ve = n;
				break;
			}
		}
	}
	assert(found);

	if (transpose)
		transpose_Matrix_d((struct Matrix_d*)xyz_ve,true);
	return ind_ve;
}
