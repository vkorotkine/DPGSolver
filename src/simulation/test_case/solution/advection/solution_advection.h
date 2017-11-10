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

#ifndef DPG__solution_advection_h__INCLUDED
#define DPG__solution_advection_h__INCLUDED
/** \file
 *  \brief Provides functions relating to the linear Advection solutions.
 */

#include "definitions_core.h"

struct Test_Case;
struct Simulation;

/// \brief Container for solution data relating to linear advection solutions.
struct Sol_Data__Advection {
	double b_adv[DMAX]; ///< The constant advection velocity vector.
};

/// \brief Set the solution function pointer members of an Euler \ref Test_Case.
void set_function_pointers_solution_advection
	(struct Test_Case* test_case,      ///< \ref Test_Case.
	 const struct Simulation*const sim ///< \ref Simulation.
	);

/** \brief Return the statically allocated \ref Sol_Data__Advection container.
 *  \return See brief. */
struct Sol_Data__Advection get_sol_data_advection
	(const struct Simulation* sim ///< \ref Simulation.
	);

/// \brief Read the required solution data into the \ref Sol_Data__Advection container.
void read_data_advection
	(const char*const input_path,              ///< Defined in \ref fopen_input.
	 struct Sol_Data__Advection*const sol_data ///< \ref Sol_Data__Advection.
	);

#endif // DPG__solution_advection_h__INCLUDED
