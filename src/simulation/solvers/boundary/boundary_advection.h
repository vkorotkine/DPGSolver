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

#ifndef DPG__boundary_advection_h__INCLUDED
#define DPG__boundary_advection_h__INCLUDED
/** \file
 *  \brief Provides containers and functions relating to boundary conditions for the linear advection equation.
 */

struct Boundary_Value_Input;
struct Boundary_Value;
struct Solver_Face;
struct Simulation;

/** \brief Version of \ref constructor_Boundary_Value_fptr computing members using the inflow (boundary) values.
 *  \return See brief. */
void constructor_Boundary_Value_advection_inflow
	(struct Boundary_Value* bv,               ///< Defined for \ref constructor_Boundary_Value_fptr.
	 const struct Boundary_Value_Input* bv_i, ///< Defined for \ref constructor_Boundary_Value_fptr.
	 const struct Solver_Face* face,          ///< Defined for \ref constructor_Boundary_Value_fptr.
	 const struct Simulation* sim             ///< Defined for \ref constructor_Boundary_Value_fptr.
	);

/** \brief Version of \ref constructor_Boundary_Value_fptr computing members using the outflow (extrapolated) values.
 *  \return See brief. */
void constructor_Boundary_Value_advection_outflow
	(struct Boundary_Value* bv,               ///< Defined for \ref constructor_Boundary_Value_fptr.
	 const struct Boundary_Value_Input* bv_i, ///< Defined for \ref constructor_Boundary_Value_fptr.
	 const struct Solver_Face* face,          ///< Defined for \ref constructor_Boundary_Value_fptr.
	 const struct Simulation* sim             ///< Defined for \ref constructor_Boundary_Value_fptr.
	);

#endif // DPG__boundary_advection_h__INCLUDED
