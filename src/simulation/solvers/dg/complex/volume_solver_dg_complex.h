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

#ifndef DPG__volume_solver_dg_complex_h__INCLUDED
#define DPG__volume_solver_dg_complex_h__INCLUDED
/** \file
 *  \brief Provides the interface for the \ref Complex_DG_Solver_Volume container and associated functions.
 *
 *  These volumes are needed by the 'D'iscontinuous 'G'alerkin solver functions for complex step linearization testing.
 */

#include <complex.h>
#include "volume_solver_dg.h"

/// \brief Container for data relating to the complex DG solver volumes.
struct Complex_DG_Solver_Volume {
	struct Solver_Volume volume; ///< The base \ref Solver_Volume.

	struct Multiarray_c* sol_coef; ///< The coefficients of the solution.
	struct Multiarray_c* rhs;      ///< The rhs terms.
};

/// \brief Constructor for a derived \ref Complex_DG_Solver_Volume.
void constructor_derived_Complex_DG_Solver_Volume
	(struct Volume* volume_ptr,   ///< Pointer to the volume.
	 const struct Simulation* sim ///< \ref Simulation.
	);

/// \brief Destructor for a derived \ref Complex_DG_Solver_Volume.
void destructor_derived_Complex_DG_Solver_Volume
	(struct Volume* volume_ptr ///< Pointer to the volume.
	);

#endif // DPG__volume_solver_dg_complex_h__INCLUDED
