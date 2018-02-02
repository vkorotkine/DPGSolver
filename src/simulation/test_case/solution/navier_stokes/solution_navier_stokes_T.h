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
 *  \brief Provides templated functions relating to the Navier-Stokes solutions.
 */

struct Test_Case_T;
struct Simulation;
struct Multiarray_T;
struct const_Multiarray_T;

/// \brief Set the solution function pointer members of an Navier-Stokes \ref Test_Case_T.
void set_function_pointers_solution_navier_stokes_T
	(struct Test_Case_T*const test_case, ///< \ref Test_Case_T.
	 const struct Simulation*const sim   ///< \ref Simulation.
	);

/// \brief Convert gradients between supported variable types (see \ref convert_variables_T for options).
void convert_variables_gradients_T
	(struct Multiarray_T*const grad,            ///< The container holding the gradient data.
	 const struct const_Multiarray_T*const sol, ///< The container holding the solution data.
	 const char type_i,                         ///< The input variable type.
	 const char type_o                          ///< The output variable type.
	);
