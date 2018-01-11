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

#ifndef DPG__definitions_geometry_h__INCLUDED
#define DPG__definitions_geometry_h__INCLUDED
/** \file
 *  \brief Provides definitions relating to the geometry.
 */

///\{ \name The supported surface parametrization options.
#define GEOM_PRM_RADIAL_PROJ 101 ///< Projection from center through affine geometry nodes to the boundary.
#define GEOM_PRM_ARC_LENGTH  102 ///< Arc length computed using numerical integration.
#define GEOM_PRM_NORMAL_PROJ 103 ///< Projection from affine geometry nodes to the boundary in the surface normal direction.
///\}

#endif // DPG__definitions_geometry_h__INCLUDED
