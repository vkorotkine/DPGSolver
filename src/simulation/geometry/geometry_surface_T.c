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
#include <string.h>
#include "gsl/gsl_math.h"

#include "macros.h"
#include "definitions_bc.h"
#include "definitions_core.h"
#include "definitions_geometry.h"
#include "definitions_mesh.h"
#include "definitions_tol.h"


#include "def_templates_geometry.h"
#include "def_templates_volume_solver.h"
#include "def_templates_matrix.h"
#include "def_templates_multiarray.h"
#include "def_templates_vector.h"
#include "def_templates_math_functions.h"
#include "def_templates_operators.h"
#include "def_templates_test_case.h"

// Static function declarations ************************************************************************************* //

/** \brief Compute the radius of the curved boundary as that of the first vertex on the boundary.
 *  \return See brief. */
static double compute_radius_from_xyz_ve
	(const int n_component_radius,                      /**< The number of components in xyz_ve to use for the
	                                                     *   computation of the radius. */
	 const struct const_Matrix_R*const xyz,             ///< The xyz coordinates in which to search.
	 const struct const_Matrix_R*const xyz_ve,          ///< The xyz coordinates of the vertices.
	 const struct Blended_Parametric_Data_T*const b_p_d ///< \ref Blended_Parametric_Data_T.
	);

// Interface functions ********************************************************************************************** //

constructor_xyz_surface_fptr_T set_constructor_xyz_surface_fptr_T
	(const char*const geom_type, const int geom_prm_type, const int domain_type)
{
	constructor_xyz_surface_fptr_T ptr = NULL;

	if (domain_type == DOM_PARAMETRIC) {
		ptr = constructor_xyz_surface_mapped_T;
	} else if (domain_type == DOM_BLENDED) {
		if (strcmp(geom_type,"n-cylinder") == 0) {
			if (geom_prm_type == GEOM_PRM_RADIAL_PROJ)
				ptr = constructor_xyz_surface_cylinder_radial_proj_T;
			else if (geom_prm_type == GEOM_PRM_ARC_LENGTH )
				ptr = constructor_xyz_surface_cylinder_arc_length_T;
			else if (geom_prm_type == GEOM_PRM_NORMAL_PROJ)
				ptr = constructor_xyz_surface_cylinder_normal_proj_T;
			else
				EXIT_ERROR("Unsupported: %d\n",geom_prm_type);
		} else {
			EXIT_ERROR("Unsupported: %s\n",geom_type);
		}
	} else {
		EXIT_ERROR("Unsupported: %d\n",domain_type);
	}
	return ptr;
}

const struct const_Matrix_T* constructor_xyz_surface_mapped_T
	(const struct Blended_Parametric_Data_T*const b_p_d)
{
	assert(b_p_d->n_type == 'c');

	const struct const_Multiarray_T*const xyz_ve_fcc = constructor_mm_NN1_Operator_const_Multiarray_T_Multiarray_R
		(b_p_d->vv0_vv_fcc,b_p_d->xyz_ve,'C','d',b_p_d->xyz_ve->order,NULL); // destructed

	const struct const_Multiarray_T*const xyz_fcc_Ma = b_p_d->constructor_xyz(0,xyz_ve_fcc,NULL,NULL); // destructed
	destructor_const_Multiarray_T(xyz_ve_fcc);

	const ptrdiff_t ext_0 = xyz_fcc_Ma->extents[0],
	                ext_1 = xyz_fcc_Ma->extents[1];
	const struct const_Matrix_T*const xyz_surf = constructor_default_const_Matrix_T(); // returned
	reinterpret_const_Multiarray_as_Matrix_T(xyz_fcc_Ma,xyz_surf,ext_0,ext_1);
	const_cast_b(&xyz_surf->owns_data,true);

	const_cast_b(&xyz_fcc_Ma->owns_data,false);
	destructor_const_Multiarray_T(xyz_fcc_Ma);

	return xyz_surf;
}

const struct const_Matrix_T* constructor_xyz_surface_cylinder_radial_proj_T
	(const struct Blended_Parametric_Data_T*const b_p_d)
{
	const struct const_Matrix_R xyz_ve_M = interpret_const_Multiarray_as_Matrix_R(b_p_d->xyz_ve);
	const struct const_Matrix_R*const xyz_ve_bX =
		constructor_mm_const_Matrix_R('N','N',1.0,b_p_d->vv0_vv_bX->op_std,&xyz_ve_M,'R'); // destructed

	const struct const_Matrix_R* xyz_ve_b = NULL;
	switch (b_p_d->n_type) {
	case 'g': // fallthrough
	case 'v':
		xyz_ve_b = xyz_ve_bX;
		break;
	case 'c': {
		xyz_ve_b = constructor_mm_const_Matrix_R('N','N',1.0,b_p_d->vv0_vv_fcc->op_std,&xyz_ve_M,'R'); // dest.
		break;
	} default:
		EXIT_ERROR("Unsupported: %c.\n",b_p_d->n_type);
		break;
	}

	const double r = compute_radius_from_xyz_ve(2,&xyz_ve_M,xyz_ve_bX,b_p_d);

	const ptrdiff_t n_n = xyz_ve_b->ext_0,
	                dim = xyz_ve_b->ext_1;
	assert(dim >= 2);

	struct Matrix_T* xyz_surf = constructor_empty_Matrix_T('R',n_n,dim); // returned
	for (int n = 0; n < n_n; ++n) {
		const double*const xyz_b = get_row_const_Matrix_R(n,xyz_ve_b);
		const double theta = atan2(xyz_b[1],xyz_b[0]);

		Type*const xyz_s = get_row_Matrix_T(n,xyz_surf);
		xyz_s[0] = r*cos(theta);
		xyz_s[1] = r*sin(theta);
		if (dim == DMAX)
			xyz_s[2] = xyz_b[2];
	}
	if (xyz_ve_b != xyz_ve_bX)
		destructor_const_Matrix_R(xyz_ve_b);
	destructor_const_Matrix_R(xyz_ve_bX);

	transpose_Matrix_T(xyz_surf,true);
	return (struct const_Matrix_T*) xyz_surf;
}

const struct const_Matrix_T* constructor_xyz_surface_cylinder_arc_length_T
	(const struct Blended_Parametric_Data_T*const b_p_d)
{
// test for arc length: Should be the same as \theta parametrization for circle.
//	1. Interpolate xyz_ve to boundary (vv0_vv_bv);
//	2. Compute parametric coordinates of vertices
//	3. Interpolate parametric coordinates to curved boundary nodes (vv0_bv_bgc)
//	4. Compute xyz curved
UNUSED(b_p_d); EXIT_ADD_SUPPORT;
}

const struct const_Matrix_T* constructor_xyz_surface_cylinder_normal_proj_T
	(const struct Blended_Parametric_Data_T*const b_p_d)
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
	(const struct const_Matrix_R*const xyz,   ///< The xyz coordinates in which to search.
	 const struct const_Matrix_R*const xyz_ve ///< The xyz coordinates of the vertices.
	);

static double compute_radius_from_xyz_ve
	(const int n_component_radius, const struct const_Matrix_R*const xyz, const struct const_Matrix_R*const xyz_ve,
	 const struct Blended_Parametric_Data_T*const b_p_d)
{
	const ptrdiff_t ind_b = find_boundary_vertex_index(xyz,xyz_ve);

	double r = 0.0;
	switch (b_p_d->domain_type) {
	case DOM_BLENDED:
		r = norm_R(n_component_radius,get_row_const_Matrix_R(ind_b,xyz),"L2");
		break;
	case DOM_PARAMETRIC: {
		const struct const_Matrix_T*const xyz_T = constructor_copy_const_Matrix_T_Matrix_R(xyz); // destructed
		const ptrdiff_t extents[] = { 1, xyz_T->ext_1, };
		const Type*const data = get_row_const_Matrix_T(ind_b,xyz_T);
		const struct const_Multiarray_T xyz_Ma = { .order = 2, .extents = extents, .layout = 'C', .data = data, };
		const struct const_Multiarray_T* xyz_surf = b_p_d->constructor_xyz(0,&xyz_Ma,NULL,NULL); // destructed
		destructor_const_Matrix_T(xyz_T);

		r = norm_R_from_T(n_component_radius,xyz_surf->data,"L2");
		destructor_const_Multiarray_T(xyz_surf);
		break;
	} default:
		EXIT_ERROR("Unsupported: %d\n",b_p_d->domain_type);
		break;
	}
	return r;
}

// Level 1 ********************************************************************************************************** //

static ptrdiff_t find_boundary_vertex_index
	(const struct const_Matrix_R*const xyz, const struct const_Matrix_R*const xyz_ve)
{
	bool found = false;
	ptrdiff_t ind_ve = -1;

	const ptrdiff_t n_n  = xyz->ext_0,
	                n_ve = xyz_ve->ext_0,
	                dim  = xyz->ext_1;
	assert(dim == xyz_ve->ext_1);

	const bool transpose = (xyz_ve->layout != 'R');
	if (transpose)
		transpose_Matrix_R((struct Matrix_R*)xyz_ve,true);

	for (int n = 0; !found && n < n_n; ++n) {
		const Real*const data_xyz = get_row_const_Matrix_R(n,xyz);
		for (int ve = 0; ve < n_ve; ++ve) {
			const Real*const data_xyz_ve = get_row_const_Matrix_R(ve,xyz_ve);
			const Real diff = norm_diff_R(dim,data_xyz_ve,data_xyz,"Inf");
			if (diff < EPS) {
				found = true;
				ind_ve = n;
				break;
			}
		}
	}
	assert(found);

	if (transpose)
		transpose_Matrix_R((struct Matrix_R*)xyz_ve,true);
	return ind_ve;
}

#include "undef_templates_geometry.h"
#include "undef_templates_volume_solver.h"
#include "undef_templates_matrix.h"
#include "undef_templates_multiarray.h"
#include "undef_templates_vector.h"
#include "undef_templates_math_functions.h"
#include "undef_templates_operators.h"
#include "undef_templates_test_case.h"
