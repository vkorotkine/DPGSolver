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

#ifndef DPG__solver_face_h__INCLUDED
#define DPG__solver_face_h__INCLUDED
/** \file
 *  \brief Provides the interface for the \ref Solver_Face container and associated functions.
 */

#include "face.h"

/// \brief Container for data relating to the solver faces.
struct Solver_Face {
	struct Face face; ///< The base \ref Face.

	/// The reference order of the face. Need not be equal to the order of the solution in the face.
	const int p_ref;
};

/** \brief Constructs the \ref Solver_Face \ref Intrusive_List.
 *  \return Standard. */
struct Intrusive_List* constructor_Solver_Faces
	(struct Simulation*const sim ///< The \ref Simulation.
	);

/// \brief Destructs the \ref Solver_Face \ref Intrusive_List.
void destructor_Solver_Faces
	(struct Intrusive_List* solver_faces ///< Standard.
	);

#endif // DPG__solver_face_h__INCLUDED
