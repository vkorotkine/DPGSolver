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
 *  \brief Provides the interface to templated functions used for surface geometry processing.
 */

#include "def_templates_geometry.h"
#include "def_templates_matrix.h"

struct Blended_Parametric_Data_T;

/** \brief Pointer to functions constructing the boundary xyz coordinates from the required members of the \ref
 *         Blended_Parametric_Data_T container for several surface parametrizations.
 *
 *  \param b_p_d \ref Blended_Parametric_Data_T.
 */
typedef const struct const_Matrix_T* (*constructor_xyz_surface_fptr_T)
	(const struct Blended_Parametric_Data_T*const b_p_d
	);


/** \brief Set the pointer to the appropriate \ref constructor_xyz_surface_fptr_T funtion based on the input geometry
 *         and surface parametrization types.
 *  \return See brief. */
constructor_xyz_surface_fptr_T set_constructor_xyz_surface_fptr_T
	(const char*const geom_type, ///< The geometry type.
	 const int geom_prm_type,    ///< The geometry parametrization type. Options: see \ref definitions_geometry.h.
	 const int domain_type       ///< \ref Simulation::domain_type.
	);

/** \brief Version of \ref constructor_xyz_surface_fptr_T using the parametric mapping.
 *  \return See brief. */
const struct const_Matrix_T* constructor_xyz_surface_mapped_T
	(const struct Blended_Parametric_Data_T*const b_p_d ///< See brief.
	);

/** \brief Version of \ref constructor_xyz_surface_fptr_T for n-cylinder surfaces (\ref GEOM_PRM_RADIAL_PROJ).
 *  \return See brief. */
const struct const_Matrix_T* constructor_xyz_surface_cylinder_radial_proj_T
	(const struct Blended_Parametric_Data_T*const b_p_d ///< See brief.
	);

/** \brief Version of \ref constructor_xyz_surface_fptr_T for n-cylinder surfaces (\ref GEOM_PRM_ARC_LENGTH).
 *  \return See brief. */
const struct const_Matrix_T* constructor_xyz_surface_cylinder_arc_length_T
	(const struct Blended_Parametric_Data_T*const b_p_d ///< See brief.
	);

/** \brief Version of \ref constructor_xyz_surface_fptr_T for n-cylinder surfaces (\ref GEOM_PRM_NORMAL_PROJ).
 *  \return See brief. */
const struct const_Matrix_T* constructor_xyz_surface_cylinder_normal_proj_T
	(const struct Blended_Parametric_Data_T*const b_p_d ///< See brief.
	);

#include "undef_templates_geometry.h"
#include "undef_templates_matrix.h"
