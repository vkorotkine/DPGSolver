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
#define DOM_STRAIGHT   1
#define DOM_CURVED     2
#define DOM_PARAMETRIC 3
///\}

///\{ \name The node coordinate tolerance of the mesh vertices
#define NODETOL_MESH 1.0e-5
///\}

#endif // DPG__definitions_mesh_h__INCLUDED
