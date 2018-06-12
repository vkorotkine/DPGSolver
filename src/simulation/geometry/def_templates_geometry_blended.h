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
 *  \brief Provides the macro definitions used for c-style templating related to the blended geometry functions.
 */

#if TYPE_RC == TYPE_REAL

///\{ \name Data types
#define Boundary_Comp_Elem_Data_T Boundary_Comp_Elem_Data_d
#define Blended_Parametric_Data_T Blended_Parametric_Data_d
///\}

///\{ \name Function names
#define constructor_xyz_blended_T                    constructor_xyz_blended
#define constructor_static_Boundary_Comp_Elem_Data_T constructor_static_Boundary_Comp_Elem_Data
#define destructor_static_Boundary_Comp_Elem_Data_T  destructor_static_Boundary_Comp_Elem_Data
#define set_Boundary_Comp_Elem_operators_T           set_Boundary_Comp_Elem_operators
#define constructor_xyz_surf_diff_T                  constructor_xyz_surf_diff
#define correct_internal_xyz_blended_T               correct_internal_xyz_blended
///\}

#define constructor_xyz_blended_ce constructor_xyz_blended_ce
#define compute_p_base_min         compute_p_base_min
#define constructor_blend_values   constructor_blend_values
#define constructor_xyz_diff_T     constructor_xyz_diff_T

#elif TYPE_RC == TYPE_COMPLEX

///\{ \name Data types
#define Boundary_Comp_Elem_Data_T Boundary_Comp_Elem_Data_c
#define Blended_Parametric_Data_T Blended_Parametric_Data_c
///\}

///\{ \name Function names
#define constructor_xyz_blended_T                    constructor_xyz_blended_c
#define constructor_static_Boundary_Comp_Elem_Data_T constructor_static_Boundary_Comp_Elem_Data_c
#define destructor_static_Boundary_Comp_Elem_Data_T  destructor_static_Boundary_Comp_Elem_Data_c
#define set_Boundary_Comp_Elem_operators_T           set_Boundary_Comp_Elem_operators_c
#define constructor_xyz_surf_diff_T                  constructor_xyz_surf_diff_c
#define correct_internal_xyz_blended_T               correct_internal_xyz_blended_c
///\}

#define constructor_xyz_blended_ce constructor_xyz_blended_ce_c
#define compute_p_base_min         compute_p_base_min_c
#define constructor_blend_values   constructor_blend_values_c
#define constructor_xyz_diff_T     constructor_xyz_diff_T_c

#endif
