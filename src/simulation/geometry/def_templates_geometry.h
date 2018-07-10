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
///\}

///\{ \name Function names
#define set_up_solver_geometry_T  set_up_solver_geometry
#define set_up_solver_geometry_p1_T set_up_solver_geometry_p1
#define compute_unit_normals_T    compute_unit_normals
#define compute_geometry_volume_T compute_geometry_volume
#define compute_NURBS_geometry_volume_T compute_NURBS_geometry_volume
#define compute_geometry_face_T   compute_geometry_face
#define compute_NURBS_geometry_face_T   compute_NURBS_geometry_face
#define constructor_xyz_s_ho_T     constructor_xyz_s_ho
#define constructor_geom_coef_ho_T constructor_geom_coef_ho
#define correct_for_exact_normals_T correct_for_exact_normals
///\}

///\{ \name Static names
#define compute_geom_coef_fptr_T compute_geom_coef_fptr_T
#define compute_normals_T compute_normals_T
#define set_fptr_geom_coef_T set_fptr_geom_coef_T
#define compute_geom_coef_p1 compute_geom_coef_p1
#define set_jacobian_permutation set_jacobian_permutation
#define compute_detJV_T compute_detJV_T
#define compute_cofactors_T compute_cofactors_T
#define constructor_xyz_fc constructor_xyz_fc
#define constructor_NURBS_xyz_fc constructor_NURBS_xyz_fc
#define constructor_xyz_fc_on_exact_boundary constructor_xyz_fc_on_exact_boundary
#define compute_unit_normals_and_det_T compute_unit_normals_and_det_T
#define compute_vol_jacobian_det_fc_T compute_vol_jacobian_det_fc_T
#define compute_geometry_volume_p1_T compute_geometry_volume_p1_T
#define compute_geometry_face_p1_T compute_geometry_face_p1_T
#define compute_geom_coef_straight_T compute_geom_coef_straight_T
#define compute_geom_coef_blended_T compute_geom_coef_blended_T
#define compute_geom_coef_parametric_T compute_geom_coef_parametric_T
#define correct_face_xyz_straight_T correct_face_xyz_straight_T
#define constructor_f_geom_coef_p1 constructor_f_geom_coef_p1
#define constructor_cc0_vgc_fgc_indices constructor_cc0_vgc_fgc_indices
#define get_operator__cv0_vgc_fis get_operator__cv0_vgc_fis
#define get_operator__vc0_fis_fgc get_operator__vc0_fis_fgc
#define get_operator__cc0_vgc_fgc get_operator__cc0_vgc_fgc
///\}

#elif TYPE_RC == TYPE_COMPLEX

///\{ \name Function pointers
#define constructor_xyz_fptr_T constructor_xyz_fptr_c
///\}

///\{ \name Function names
#define set_up_solver_geometry_T  set_up_solver_geometry_c
#define set_up_solver_geometry_p1_T set_up_solver_geometry_p1_c
#define compute_unit_normals_T    compute_unit_normals_c
#define compute_geometry_volume_T compute_geometry_volume_c
#define compute_NURBS_geometry_volume_T compute_NURBS_geometry_volume_c
#define compute_geometry_face_T   compute_geometry_face_c
#define compute_NURBS_geometry_face_T   compute_NURBS_geometry_face_c
#define constructor_xyz_s_ho_T     constructor_xyz_s_ho_c
#define constructor_geom_coef_ho_T constructor_geom_coef_ho_c
#define correct_for_exact_normals_T correct_for_exact_normals_c
///\}

///\{ \name Static names
#define compute_geom_coef_fptr_T compute_geom_coef_fptr_T_c
#define compute_normals_T compute_normals_T_c
#define set_fptr_geom_coef_T set_fptr_geom_coef_T_c
#define compute_geom_coef_p1 compute_geom_coef_p1_c
#define set_jacobian_permutation set_jacobian_permutation_c
#define compute_detJV_T compute_detJV_T_c
#define compute_cofactors_T compute_cofactors_T_c
#define constructor_xyz_fc constructor_xyz_fc_c
#define constructor_NURBS_xyz_fc constructor_NURBS_xyz_fc_c
#define constructor_xyz_fc_on_exact_boundary constructor_xyz_fc_on_exact_boundary_c
#define compute_unit_normals_and_det_T compute_unit_normals_and_det_T_c
#define compute_vol_jacobian_det_fc_T compute_vol_jacobian_det_fc_T_c
#define compute_geometry_volume_p1_T compute_geometry_volume_p1_T_c
#define compute_geometry_face_p1_T compute_geometry_face_p1_T_c
#define compute_geom_coef_straight_T compute_geom_coef_straight_T_c
#define compute_geom_coef_blended_T compute_geom_coef_blended_T_c
#define compute_geom_coef_parametric_T compute_geom_coef_parametric_T_c
#define correct_face_xyz_straight_T correct_face_xyz_straight_T_c
#define constructor_f_geom_coef_p1 constructor_f_geom_coef_p1_c
#define constructor_cc0_vgc_fgc_indices constructor_cc0_vgc_fgc_indices_c
#define get_operator__cv0_vgc_fis get_operator__cv0_vgc_fis_c
#define get_operator__vc0_fis_fgc get_operator__vc0_fis_fgc_c
#define get_operator__cc0_vgc_fgc get_operator__cc0_vgc_fgc_c
///\}

#endif
