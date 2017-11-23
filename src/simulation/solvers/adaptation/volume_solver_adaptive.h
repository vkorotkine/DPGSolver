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

#ifndef DPG__volume_solver_adaptive_h__INCLUDED
#define DPG__volume_solver_adaptive_h__INCLUDED
/** \file
 *  \brief Provides the interface for the \ref Adaptive_Solver_Volume container and associated functions.
 */

#include "volume_solver.h"

/// \brief Container for data relating to the adaptive solver volumes.
struct Adaptive_Solver_Volume {
	struct Solver_Volume volume; ///< The base \ref Solver_Volume.

	int adapt_type; ///< The type of adaptation to use. Options: see \ref definitions_adaption.h.

	int p_ref_prev; ///< The previous value of \ref Solver_Volume::p_ref.

	bool updated; ///< Flag for whether the volume was updated.
};

/// \brief Constructor for a derived \ref Adaptive_Solver_Volume.
void constructor_derived_Adaptive_Solver_Volume
	(struct Volume* volume_ptr,   ///< Pointer to the volume.
	 const struct Simulation* sim ///< \ref Simulation.
	);

/// \brief Destructor for a derived \ref Adaptive_Solver_Volume.
void destructor_derived_Adaptive_Solver_Volume
	(struct Volume* volume_ptr ///< Pointer to the volume.
	);

#endif // DPG__volume_solver_adaptive_h__INCLUDED
