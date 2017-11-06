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

#ifndef DPG__face_solver_dg_complex_h__INCLUDED
#define DPG__face_solver_dg_complex_h__INCLUDED
/** \file
 *  \brief Provides the interface for the \ref Complex_DG_Solver_Face container and associated functions.
 *
 *  These faces are needed by the 'D'iscontinuous 'G'alerkin solver functions for complex step linearization testing.
 */

#include <complex.h>
#include "face_solver_dg.h"
#include "complex_boundary.h"

/// \brief Container for data relating to the complex DG solver faces.
struct Complex_DG_Solver_Face {
	struct DG_Solver_Face face; ///< The base \ref DG_Solver_Face.

	/// Complex version of \ref constructor_Boundary_Value_fptr.
	constructor_Boundary_Value_c_fptr constructor_Boundary_Value_c_fcl;
};

/// \brief Constructor for a derived \ref Complex_DG_Solver_Face.
void constructor_derived_Complex_DG_Solver_Face
	(struct Face* face_ptr,       ///< Pointer to the face.
	 const struct Simulation* sim ///< \ref Simulation.
	);

/// \brief Destructor for a derived \ref Complex_DG_Solver_Face.
void destructor_derived_Complex_DG_Solver_Face
	(struct Face* face_ptr ///< Pointer to the face.
	);

#endif // DPG__face_solver_dg_complex_h__INCLUDED
