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

///\{ \name Function pointers.
#define constructor_xyz_surface_fptr_T constructor_xyz_surface_fptr
///\}

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

#define set_constructor_xyz_surface_fptr_T             set_constructor_xyz_surface_fptr_d
#define constructor_xyz_surface_cylinder_radial_proj_T constructor_xyz_surface_cylinder_radial_proj_d
#define constructor_xyz_surface_cylinder_arc_length_T  constructor_xyz_surface_cylinder_arc_length_d
#define constructor_xyz_surface_cylinder_normal_proj_T constructor_xyz_surface_cylinder_normal_proj_d
///\}

#elif TYPE_RC == TYPE_COMPLEX

///\{ \name Function pointers.
#define constructor_xyz_surface_fptr_T constructor_xyz_surface_fptr_c
///\}

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

#define set_constructor_xyz_surface_fptr_T             set_constructor_xyz_surface_fptr_c
#define constructor_xyz_surface_cylinder_radial_proj_T constructor_xyz_surface_cylinder_radial_proj_c
#define constructor_xyz_surface_cylinder_arc_length_T  constructor_xyz_surface_cylinder_arc_length_c
#define constructor_xyz_surface_cylinder_normal_proj_T constructor_xyz_surface_cylinder_normal_proj_c
///\}

#endif
