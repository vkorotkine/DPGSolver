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
 *  \brief Provides the interface to templated functions used to solve for the solution.
 */

#include <stddef.h>

struct Simulation;

/** \brief Constructor for a \ref Solver_Storage_Implicit container.
 *  \return See brief. */
struct Solver_Storage_Implicit* constructor_Solver_Storage_Implicit_T
	(const struct Simulation* sim ///< \ref Simulation.
	);

/** \brief Compute the number of 'd'egrees 'o'f 'f'reedom.
 *  \return See brief. */
ptrdiff_t compute_dof_T
	(const struct Simulation* sim ///< \ref Simulation.
	);

/// \brief Update \ref Solver_Volume_T::ind_dof and \ref Solver_Face_T::ind_dof.
void update_ind_dof_T
	(const struct Simulation* sim ///< \ref Simulation.
	);
