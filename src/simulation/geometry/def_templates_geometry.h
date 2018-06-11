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
 *  \brief Provides the macro definitions used for c-style templating related to the geometry functions.
 */

#include "def_templates_geometry_blended.h"
#include "def_templates_geometry_normals.h"
#include "def_templates_geometry_parametric.h"
#include "def_templates_geometry_surface.h"

#if TYPE_RC == TYPE_REAL

///\{ \name Function pointers
#define constructor_xyz_fptr_T constructor_xyz_fptr

#define compute_geom_coef_fptr_T compute_geom_coef_fptr
///\}

///\{ \name Function names
#define set_up_solver_geometry_T  set_up_solver_geometry
#define set_up_solver_geometry_p1_T set_up_solver_geometry_p1
#define compute_unit_normals_T    compute_unit_normals
#define compute_geometry_volume_T compute_geometry_volume
#define compute_geometry_face_T   compute_geometry_face
#define constructor_xyz_s_ho_T     constructor_xyz_s_ho
#define constructor_geom_coef_ho_T constructor_geom_coef_ho
#define correct_for_exact_normals_T correct_for_exact_normals
///\}

#elif TYPE_RC == TYPE_COMPLEX

///\{ \name Function pointers
#define constructor_xyz_fptr_T constructor_xyz_fptr_c

#define compute_geom_coef_fptr_T compute_geom_coef_fptr_c
///\}

///\{ \name Function names
#define set_up_solver_geometry_T  set_up_solver_geometry_c
#define set_up_solver_geometry_p1_T set_up_solver_geometry_p1_c
#define compute_unit_normals_T    compute_unit_normals_c
#define compute_geometry_volume_T compute_geometry_volume_c
#define compute_geometry_face_T   compute_geometry_face_c
#define constructor_xyz_s_ho_T     constructor_xyz_s_ho_c
#define constructor_geom_coef_ho_T constructor_geom_coef_ho_c
#define correct_for_exact_normals_T correct_for_exact_normals_c
///\}

#endif
