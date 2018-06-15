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
 *  \brief Provides the interface to templated functions used for exact normals.
 */

#include "def_templates_geometry.h"
#include "def_templates_face_solver.h"
#include "def_templates_multiarray.h"

struct Solver_Face_T;
struct Simulation;

/** \brief Function pointer to exact boundary normal computing functions.
 *  \return Standard.
 *
 *  \param s_face \ref Solver_Face_T.
 */
typedef void (*correct_for_exact_normal_fptr_T)
	(struct Solver_Face_T*const s_face
	);

/** \brief Return the appropriate \ref correct_for_exact_normal_fptr_T pointer based on the test case.
 *  \return See brief. */
correct_for_exact_normal_fptr_T set_correct_for_exact_normal_fptr_T
	(const struct Simulation*const sim ///< Standard.
	);

/// \brief Version of \ref correct_for_exact_normal_fptr_T for cylinder geometry.
void correct_for_exact_normal_cylinder_T
	(struct Solver_Face_T*const s_face ///< See brief.
	);

/// \brief Version of \ref correct_for_exact_normal_fptr_T for Gaussian bump geometry.
void correct_for_exact_normal_gaussian_bump_T
	(struct Solver_Face_T*const s_face ///< See brief.
	);

#include "undef_templates_geometry.h"
#include "undef_templates_face_solver.h"
#include "undef_templates_multiarray.h"
