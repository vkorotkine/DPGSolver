// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#ifndef DPG__mesh_geometry_cylinder_hollow_section_h__INCLUDED
#define DPG__mesh_geometry_cylinder_hollow_section_h__INCLUDED
/**	\file
 *	Provides functions relating to hollow cylindrical geometry.
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
