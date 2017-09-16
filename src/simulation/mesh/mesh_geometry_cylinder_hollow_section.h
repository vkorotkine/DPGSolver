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

#ifndef DPG__mesh_geometry_cylinder_hollow_section_h__INCLUDED
#define DPG__mesh_geometry_cylinder_hollow_section_h__INCLUDED
/**	\file
 *	\brief Provides functions relating to hollow cylindrical geometry.
 */

struct const_Vector_i;
struct Matrix_d;

/// \brief Snaps vertices to the cylinder hollow section geometry.
void mesh_snap_to_cylinder__hollow_section
	(const char*const input_path,                 ///< Defined in \ref mesh_snap_to_boundary_fptr.
	 const struct const_Vector_i*const ve_curved, ///< Defined in \ref mesh_snap_to_boundary_fptr.
	 const struct Matrix_d*const nodes            ///< Defined in \ref mesh_snap_to_boundary_fptr.
	);

#endif // DPG__mesh_geometry_cylinder_hollow_section_h__INCLUDED
