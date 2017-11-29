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

#ifndef DPG__test_complex_boundary_euler_h__INCLUDED
#define DPG__test_complex_boundary_euler_h__INCLUDED
/** \file
 *  \brief Provides `complex` versions of containers and functions defined in \ref boundary_euler.h.
 */

struct Boundary_Value_Input_c;
struct Boundary_Value_c;
struct Solver_Face;
struct Simulation;

/** \brief `complex` version of \ref constructor_Boundary_Value_euler_riemann.
 *  \return See brief. */
void constructor_Boundary_Value_c_euler_riemann
	(struct Boundary_Value_c* bv,               ///< See brief.
	 const struct Boundary_Value_Input_c* bv_i, ///< See brief.
	 const struct Solver_Face* face,            ///< See brief.
	 const struct Simulation* sim               ///< See brief.
	);

/** \brief `complex` version of \ref constructor_Boundary_Value_euler_slipwall.
 *  \return See brief. */
void constructor_Boundary_Value_c_euler_slipwall
	(struct Boundary_Value_c* bv,               ///< See brief.
	 const struct Boundary_Value_Input_c* bv_i, ///< See brief.
	 const struct Solver_Face* face,            ///< See brief.
	 const struct Simulation* sim               ///< See brief.
	);

#endif // DPG__test_complex_boundary_euler_h__INCLUDED
