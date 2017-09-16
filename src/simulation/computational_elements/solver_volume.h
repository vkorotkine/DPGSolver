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

#ifndef DPG__solver_volume_h__INCLUDED
#define DPG__solver_volume_h__INCLUDED
/** \file
 *  \brief Provides the interface for the \ref Solver_Volume container and associated functions.
 */

#include "volume.h"

/// \brief Container for data relating to the solver volumes.
struct Solver_Volume {
	struct Volume volume; ///< The base \ref Volume.

	int p; ///< The order of the solution.

	/// The coefficients of the solution in the \ref Simulation::basis_sol.
	struct Multiarray_d* sol_coef;

	/// The coefficients of the solution gradient in the \ref Simulation::basis_sol.
	struct Multiarray_d* grad_coef;

	/** The metric terms (cofactors of the geometry Jacobian) used for transformation of integrals between physical
	 *  and computational space stored at the (v)olume (g)eometry nodes. */
	const struct const_Multiarray_d*const metrics_vg;

	/** The metric terms (cofactors of the geometry Jacobian) used for transformation of integrals between physical
	 *  and computational space stored at the (v)olume (c)ubature nodes. */
	const struct const_Multiarray_d*const metrics_vc;

	/// The determinate of the geometry Jacobian evaluated at the volume cubature nodes.
	const struct const_Multiarray_d*const jacobian_det_vc;
};

/** \brief Constructs the \ref Solver_Volume \ref Intrusive_List.
 *  \return Standard. */
struct Intrusive_List* constructor_Solver_Volumes
	(struct Simulation*const sim ///< The \ref Simulation.
	);

/// \brief Destructs the \ref Solver_Volume \ref Intrusive_List.
void destructor_Solver_Volumes
	(struct Intrusive_List* solver_volumes ///< Standard.
	);

#endif // DPG__solver_volume_h__INCLUDED
