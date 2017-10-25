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

#ifndef DPG__compute_error_h__INCLUDED
#define DPG__compute_error_h__INCLUDED
/** \file
 *  \brief Provides the interface to functions used for error computation and output.
 */

struct Error_CE;
struct Simulation;
struct Solver_Volume;

#include "compute_error_euler.h"

/** \brief Function pointer to error computing functions (which construct \ref Error_CE containers).
 *  \param sim \ref Simulation.
 */
typedef struct Error_CE* (*constructor_Error_CE_fptr)
	(const struct Simulation* sim
	);

/// Container holding information relating to the computational elements.
struct Error_CE {
	const double domain_volume; ///< The volume of the domain.

	const struct const_Vector_d* sol_L2; ///< The solution measured in L2.
};

// Interface functions ********************************************************************************************** //

/// \brief Output the error of the solution.
void output_error
	(const struct Simulation* sim ///< \ref Simulation.
	);

/** \brief Compute the volume of the input \ref Solver_Volume.
 *  \return See brief. */
double compute_volume
	(const struct Solver_Volume* s_vol ///< \ref Solver_Volume.
	);

#endif // DPG__compute_error_h__INCLUDED
