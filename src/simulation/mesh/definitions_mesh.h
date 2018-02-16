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

#ifndef DPG__definitions_mesh_h__INCLUDED
#define DPG__definitions_mesh_h__INCLUDED
/**	\file
 *	\brief Provides the definitions relating to the \ref Mesh container.
 */

///\{ \name The supported domain types.
#define DOM_STRAIGHT   1 ///< Affine mesh.

/// Affine mesh excluding the boundary elements which have the curved boundary representation blended into the volume.
#define DOM_BLENDED    2

#define DOM_PARAMETRIC 3 ///< Parametric mesh. Generally non-affine.
///\}

///\{ \name The node coordinate tolerance of the mesh vertices
#define NODETOL_MESH 1.0e-5
///\}


///\{ \name "Boundary conditions" which result in no boundary (i.e. the boundaries connect to each other periodically).
#define PERIODIC_XL 51
#define PERIODIC_XR 52
#define PERIODIC_YL 53
#define PERIODIC_YR 54
#define PERIODIC_ZL 55
#define PERIODIC_ZR 56

#define PERIODIC_XL_REFLECTED_Y 61
#define PERIODIC_XR_REFLECTED_Y 62
///\}

#endif // DPG__definitions_mesh_h__INCLUDED
