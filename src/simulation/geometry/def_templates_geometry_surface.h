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
 *  \brief Provides the macro definitions used for c-style templating related to the surface geometry functions.
 */

#if TYPE_RC == TYPE_REAL

///\{ \name Function pointers.
#define constructor_xyz_surface_fptr_T constructor_xyz_surface_fptr
///\}

///\{ \name Function names
#define set_constructor_xyz_surface_fptr_T             set_constructor_xyz_surface_fptr_d
#define constructor_xyz_surface_mapped_T               constructor_xyz_surface_mapped_d
#define constructor_xyz_surface_cylinder_radial_proj_T constructor_xyz_surface_cylinder_radial_proj_d
#define constructor_xyz_surface_cylinder_arc_length_T  constructor_xyz_surface_cylinder_arc_length_d
#define constructor_xyz_surface_cylinder_normal_proj_T constructor_xyz_surface_cylinder_normal_proj_d
///\}

#elif TYPE_RC == TYPE_COMPLEX

///\{ \name Function pointers.
#define constructor_xyz_surface_fptr_T constructor_xyz_surface_fptr_c
///\}

///\{ \name Function names
#define set_constructor_xyz_surface_fptr_T             set_constructor_xyz_surface_fptr_c
#define constructor_xyz_surface_mapped_T               constructor_xyz_surface_mapped_c
#define constructor_xyz_surface_cylinder_radial_proj_T constructor_xyz_surface_cylinder_radial_proj_c
#define constructor_xyz_surface_cylinder_arc_length_T  constructor_xyz_surface_cylinder_arc_length_c
#define constructor_xyz_surface_cylinder_normal_proj_T constructor_xyz_surface_cylinder_normal_proj_c
///\}

#endif
