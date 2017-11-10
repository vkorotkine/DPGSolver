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

#ifndef DPG__volume_solver_dpg_h__INCLUDED
#define DPG__volume_solver_dpg_h__INCLUDED
/** \file
 *  \brief Provides the interface for the \ref DPG_Solver_Volume container and associated functions.
 *
 *  These volumes are needed by the 'D'iscontinuous 'P'etrov 'G'alerkin solver functions.
 */

#include "volume_solver.h"

/// \brief Container for data relating to the DPG solver volumes.
struct DPG_Solver_Volume {
	struct Solver_Volume volume; ///< The base \ref Solver_Volume.

	const struct const_Matrix_d* norm_op_H0; ///< The H0 (L2) norm operator.
};

/// \brief Constructor for a derived \ref DPG_Solver_Volume.
void constructor_derived_DPG_Solver_Volume
	(struct Volume* volume_ptr,   ///< Pointer to the volume.
	 const struct Simulation* sim ///< \ref Simulation.
	);

/// \brief Destructor for a derived \ref DPG_Solver_Volume.
void destructor_derived_DPG_Solver_Volume
	(struct Volume* volume_ptr ///< Pointer to the volume.
	);

#endif // DPG__volume_solver_dpg_h__INCLUDED
