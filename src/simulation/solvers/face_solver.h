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

#ifndef DPG__face_solver_h__INCLUDED
#define DPG__face_solver_h__INCLUDED
/** \file
 *  \brief Provides the interface for the real \ref Solver_Face_T container and associated functions.
 */

#include "face.h"
#include "boundary.h"

#include "def_templates_type_d.h"
#include "def_templates_multiarray.h"
#include "def_templates_vector.h"
#include "def_templates_boundary.h"
#include "def_templates_face_solver.h"
#include "face_solver_T.h"
#include "undef_templates_type.h"
#include "undef_templates_multiarray.h"
#include "undef_templates_vector.h"
#include "undef_templates_boundary.h"
#include "undef_templates_face_solver.h"

struct Solver_Face_c;

/// \brief Copy members from a real to a complex \ref Solver_Face_T.
void copy_members_r_to_c_Solver_Face
	(struct Solver_Face_c*const s_face,       ///< The complex \ref Solver_Face_T.
	 const struct Solver_Face*const s_face_r, ///< The real \ref Solver_Face_T.
	 const struct Simulation*const sim        ///< \ref Simulation.
	);

#endif // DPG__face_solver_h__INCLUDED
