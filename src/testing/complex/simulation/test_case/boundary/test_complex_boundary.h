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

#ifndef DPG__test_complex_boundary_h__INCLUDED
#define DPG__test_complex_boundary_h__INCLUDED
/** \file
 *  \brief Provides `complex` versions of containers and functions defined in \ref boundary.h.
 */

#include <stdbool.h>
#include "boundary.h"
#include "complex_boundary.h"

struct Boundary_Value_Input_c;
struct Boundary_Value_c;
struct Solver_Face;
struct Simulation;

/** \brief `complex` version of \ref constructor_Boundary_Value_Input_face_fptr.
 *  \return Standard.
 *
 *  \param bv_i   See brief.
 *  \param s_face See brief.
 *  \param sim    See brief.
 */
typedef void (*constructor_Boundary_Value_Input_c_face_fptr)
	(struct Boundary_Value_Input_c* bv_i,
	 const struct Solver_Face* s_face,
	 const struct Simulation* sim
	);

/// \brief Derived `complex` version of \ref Boundary_Value_Input.
struct Boundary_Value_Input_c {
	struct Boundary_Value_Input bv_i; ///< Base \ref Boundary_Value_Input.

	const struct const_Multiarray_c* s; ///< See brief.
	const struct const_Multiarray_c* g; ///< See brief.
};

/// \brief `complex` version of \ref Boundary_Value.
struct Boundary_Value_c {
	const struct const_Multiarray_c* s; ///< See brief.
	const struct const_Multiarray_c* g; ///< See brief.
};

// Interface functions ********************************************************************************************** //

/** \brief `complex` version of \ref constructor_Boundary_Value_Input_face_s_fcl_interp.
 *  \return See brief. */
void constructor_Boundary_Value_Input_c_face_s_fcl_interp
	(struct Boundary_Value_Input_c* bv_i, ///< See brief.
	 const struct Solver_Face* face,      ///< See brief.
	 const struct Simulation* sim         ///< See brief.
	);

/// \brief Destructor for a \ref Boundary_Value_Input_c container.
void destructor_Boundary_Value_Input_c
	(struct Boundary_Value_Input_c* bv_i ///< Standard.
	);

/** \brief `complex` version of \ref constructor_Boundary_Value_s_fcl_interp.
 *  \return See brief. */
void constructor_Boundary_Value_c_s_fcl_interp
	(struct Boundary_Value_c* bv,               ///< See brief.
	 const struct Boundary_Value_Input_c* bv_i, ///< See brief.
	 const struct Solver_Face* face,            ///< See brief.
	 const struct Simulation* sim               ///< See brief.
	);

/// \brief Destructor for a \ref Boundary_Value_c container.
void destructor_Boundary_Value_c
	(struct Boundary_Value_c* bv ///< Standard.
	);

#endif // DPG__test_complex_boundary_h__INCLUDED
