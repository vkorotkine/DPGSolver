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

#ifndef DPG__mesh_periodic_h__INCLUDED
#define DPG__mesh_periodic_h__INCLUDED
/**	\file
 *	\brief Provides the interface to mesh periodic containers and functions.
 */

struct Mesh_Data;
struct Conn_info;

/** \brief Correct the face vertex correspondence if the mesh is periodic.
 *
 *	The correction modifies the vertices of the "slave" faces to those of the "master" faces. The correspondence is
 *	established by checking the equivalence of the mesh node coordinates in the appropriate directions.
 */
void correct_f_ve_for_periodic
	(const struct Mesh_Data*const mesh_data, ///< The \ref Mesh_Data.
	 struct Conn_info*const conn_info        ///< The \ref Conn_info.
	);

#endif // DPG__mesh_periodic_h__INCLUDED
