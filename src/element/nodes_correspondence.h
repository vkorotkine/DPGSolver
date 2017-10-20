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

#ifndef DPG__nodes_correspondence_h__INCLUDED
#define DPG__nodes_correspondence_h__INCLUDED
/** \file
 *  \brief Provides the interface to functions computing the index correspondence between face nodes as seen from
 *         neighbouring volumes on either side.
 */

///\{ \name Number of possible permutations for two neighbouring faces.
#define POINT_N_PERM 1
#define LINE_N_PERM  2
#define TRI_N_PERM   6
#define QUAD_N_PERM  8
///\}

/** \brief Constructor for \ref const_Multiarray_Vector_i\* container holding the correspondence between node
 *         permutations on the same face as seen from the two neighbouring volumes.
 *  \return See brief.
 *
 *  The ordering is such that when you re-order the nodes using the face correspondence indices, you get back the
 *  ordering on the reference element. See the data file for \ref constructor_Nodes_FC_Data for visualizations.
 */
const struct const_Multiarray_Vector_i* constructor_nodes_face_corr
	(const int d,         ///< Defined in \ref constructor_Nodes_fptr.
	 const int p,         ///< Defined in \ref constructor_Nodes_fptr.
	 const int node_type, ///< Defined in \ref constructor_Nodes_fptr.
	 const int s_type     ///< \ref Element::s_type.
	);

#endif // DPG__nodes_correspondence_h__INCLUDED
