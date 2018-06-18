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

#ifndef DPG__test_support_solve_h__INCLUDED
#define DPG__test_support_solve_h__INCLUDED
/** \file
 *  \brief Provides supporting functions for testing of the various solvers.
 */

#include <stdbool.h>

struct Face;
struct Volume;
struct Simulation;
struct Intrusive_List;

/// \brief Perturb the initial solution parameters for the supported methods.
void perturb_solution
	(const struct Simulation*const sim ///< Standard.
		);

/** \brief Constructor for the list of \ref Volume\*s including only the current volume.
 *  \return Standard. */
struct Intrusive_List* constructor_Volumes_local_centre_only
	(const struct Volume*const vol ///< The centre \ref Volume.
		);

/** \brief Constructor for the list of \ref Face\*s adjacent to the current volume.
 *  \return Standard. */
struct Intrusive_List* constructor_Faces_local_neigh_only
	(const struct Volume*const vol,    ///< The centre \ref Volume.
	 const struct Simulation*const sim ///< Standard.
		);

/** \brief Check whether the \ref Face neighbours the current \ref Volume.
 *  \return `true` if neighbouring; `false` otherwise. */
bool is_face_neighbour
	(const struct Face*const face,      ///< The face under investigation.
	 const struct Volume*const vol_curr ///< The current volume.
		);

#endif // DPG__test_support_solve_h__INCLUDED
