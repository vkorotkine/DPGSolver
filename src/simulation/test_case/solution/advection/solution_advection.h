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
 *  \brief Provides real functions relating to the linear Advection solutions.
 */

#include "def_templates_type_d.h"
#include "def_templates_solution_advection.h"
#include "def_templates_test_case.h"
#include "solution_advection_T.h"
#include "undef_templates_type.h"
#include "undef_templates_solution_advection.h"
#include "undef_templates_test_case.h"


#include "definitions_core.h"
struct Simulation;

/// \brief Container for solution data relating to linear advection solutions.
struct Sol_Data__Advection {
	double b_adv[DMAX]; ///< The constant advection velocity vector.
};

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
