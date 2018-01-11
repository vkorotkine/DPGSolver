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

#ifndef DPG__geometry_blended_h__INCLUDED
#define DPG__geometry_blended_h__INCLUDED
/** \file
 *  \brief Provides the interface to real functions used for blended geometry processing.
 */

#include "def_templates_type_d.h"
#include "def_templates_geometry_blended.h"
#include "def_templates_volume_solver.h"
#include "def_templates_multiarray.h"
#include "geometry_blended_T.h"
#include "undef_templates_type.h"
#include "undef_templates_geometry_blended.h"
#include "undef_templates_volume_solver.h"
#include "undef_templates_multiarray.h"

struct Blended_Parametric_Data;

/** \brief Pointer to functions constructing the boundary xyz coordinates from the required members of the \ref
 *         Blended_Parametric_Data container for several surface parametrizations.
 *
 *  \param b_p_d \ref Blended_Parametric_Data.
 */
typedef const struct const_Matrix_d* (*constructor_xyz_surface_fptr)
	(const struct Blended_Parametric_Data*const b_p_d
	);

/** \brief Container for data required to compute curved surface geometry values for various methods of surface
 *         parametrization.
 */
struct Blended_Parametric_Data {
	const int geom_parametrization; ///< \ref Test_Case_T::geom_parametrization.

	constructor_xyz_surface_fptr constructor_xyz_surface; ///< \ref Solver_Volume_T::constructor_xyz_surface.

	const struct const_Multiarray_d* xyz_ve; ///< \ref Volume::xyz_ve.

	const struct Operator* vv0_vv_bX; ///< See notation in \ref element_operators.h.
	const struct Operator* vv0_vv_bv; ///< See notation in \ref element_operators.h.
	const struct Operator* vv0_bv_bX; ///< See notation in \ref element_operators.h.

	const struct const_Vector_d* normal; ///< Outward pointing unit normal vector.
};

/** \brief Set the pointer to the appropriate \ref constructor_xyz_surface_fptr funtion based on the input geometry and
 *         surface parametrization types.
 *  \return See brief. */
constructor_xyz_surface_fptr set_constructor_xyz_surface_fptr
	(const char*const geom_type, ///< The geometry type.
	 const int geom_prm_type     ///< The geometry parametrization type. Options: see \ref definitions_geometry.h.
	);

// Geometry-specific functions ************************************************************************************** //

/** \brief Version of \ref constructor_xyz_surface_fptr for n-cylinder surfaces (\ref GEOM_PRM_RADIAL_PROJ).
 *  \return See brief. */
const struct const_Matrix_d* constructor_xyz_surface_cylinder_radial_proj
	(const struct Blended_Parametric_Data*const b_p_d ///< See brief.
	);

/** \brief Version of \ref constructor_xyz_surface_fptr for n-cylinder surfaces (\ref GEOM_PRM_ARC_LENGTH).
 *  \return See brief. */
const struct const_Matrix_d* constructor_xyz_surface_cylinder_arc_length
	(const struct Blended_Parametric_Data*const b_p_d ///< See brief.
	);

/** \brief Version of \ref constructor_xyz_surface_fptr for n-cylinder surfaces (\ref GEOM_PRM_NORMAL_PROJ).
 *  \return See brief. */
const struct const_Matrix_d* constructor_xyz_surface_cylinder_normal_proj
	(const struct Blended_Parametric_Data*const b_p_d ///< See brief.
	);

#endif // DPG__geometry_blended_h__INCLUDED
