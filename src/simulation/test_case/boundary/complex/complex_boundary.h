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

#ifndef DPG__complex_boundary_h__INCLUDED
#define DPG__complex_boundary_h__INCLUDED
/** \file
 *  \brief Provides the **minimal** interface for containers and functions relating to complex boundary conditions.
 */

struct Boundary_Value_Input_c;
struct Boundary_Value_c;
struct Solver_Face;
struct Simulation;

/** \brief `complex` version of \ref constructor_Boundary_Value_fptr.
 *  \return Standard.
 *
 *  \param bv     See brief.
 *  \param bv_i   See brief.
 *  \param s_face See brief.
 *  \param sim    See brief.
 */
typedef void (*constructor_Boundary_Value_c_fptr)
	(struct Boundary_Value_c* bv,
	 const struct Boundary_Value_Input_c* bv_i,
	 const struct Solver_Face* s_face,
	 const struct Simulation* sim
	);

#endif // DPG__complex_boundary_h__INCLUDED
