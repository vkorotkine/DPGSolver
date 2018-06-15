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

#ifndef DPG__face_solver_dpg_h__INCLUDED
#define DPG__face_solver_dpg_h__INCLUDED
/** \file
 *  \brief Provides the interface for the real \ref DPG_Solver_Face_T container and associated functions.
 */

#include "face_solver.h"

#include "def_templates_type_d.h"
#include "face_solver_dpg_T.h"
#include "undef_templates_type.h"

#include "def_templates_type_dc.h"
#include "face_solver_dpg_T.h"
#include "undef_templates_type.h"

struct DPG_Solver_Face_c;

/// \brief Copy members from a real to a complex \ref DPG_Solver_Face_T.
void copy_members_r_to_c_DPG_Solver_Face
	(struct DPG_Solver_Face_c*const dpg_s_face,       ///< The complex \ref DPG_Solver_Face_T.
	 const struct DPG_Solver_Face*const dpg_s_face_r, ///< The real \ref DPG_Solver_Face_T.
	 const struct Simulation*const sim                ///< \ref Simulation.
	);

#endif // DPG__face_solver_dpg_h__INCLUDED
