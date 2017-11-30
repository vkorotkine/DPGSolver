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

#ifndef DPG__solution_euler_h__INCLUDED
#define DPG__solution_euler_h__INCLUDED
/** \file
 *  \brief Provides functions relating to the Euler solutions.
 */

struct Test_Case;
struct Simulation;

struct Multiarray_d;
struct const_Multiarray_d;

/// \brief Set the solution function pointer members of an Euler \ref Test_Case.
void set_function_pointers_solution_euler
	(struct Test_Case* test_case,      ///< \ref Test_Case.
	 const struct Simulation*const sim ///< \ref Simulation.
	);

/// \brief `double` version of \ref convert_variables_T.
void convert_variables
	(struct Multiarray_d* vars, ///< See brief.
	 const char type_i,         ///< See brief.
	 const char type_o          ///< See brief.
	);

/// \brief Compute the entropy measure: s = p/pow(rho,GAMMA).
void compute_entropy
	(struct Multiarray_d* s,                ///< The container to hold the entropy data.
	 const struct const_Multiarray_d* vars, ///< The container of Euler variables.
	 const char var_type                    ///< The type of the variables.
	);

/// \brief Compute the mach number: mach = V/c.
void compute_mach
	(struct Multiarray_d* mach,             ///< The container to hold the mach data.
	 const struct const_Multiarray_d* vars, ///< The container of Euler variables.
	 const char var_type                    ///< The type of the variables.
	);

#endif // DPG__solution_euler_h__INCLUDED
