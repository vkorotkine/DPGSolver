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
 *  \brief Provides templated containers and functions relating to boundary conditions for the
 *         linear advection equation.
 */

struct Boundary_Value_Input_T;
struct Boundary_Value_T;
struct Solver_Face_T;
struct Simulation;

/// \brief Version of \ref constructor_Boundary_Value_fptr computing members using the inflow (boundary) values.
void constructor_Boundary_Value_T_advection_inflow
	(struct Boundary_Value_T* bv,               ///< See brief.
	 const struct Boundary_Value_Input_T* bv_i, ///< See brief.
	 const struct Solver_Face_T* face,          ///< See brief.
	 const struct Simulation* sim               ///< See brief.
	);

/// \brief Version of \ref constructor_Boundary_Value_fptr computing members using the outflow (extrapolated) values.
void constructor_Boundary_Value_T_advection_outflow
	(struct Boundary_Value_T* bv,               ///< See brief.
	 const struct Boundary_Value_Input_T* bv_i, ///< See brief.
	 const struct Solver_Face_T* face,          ///< See brief.
	 const struct Simulation* sim               ///< See brief.
	);
