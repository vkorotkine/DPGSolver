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
 *  \brief Provides the interface for the templated \ref DG_Solver_Volume_T container and associated functions.
 *
 *  These volumes are needed by the 'D'iscontinuous 'G'alerkin solver functions.
 */

#include "definitions_core.h"

/// \brief Container for data relating to the DG solver volumes.
struct DG_Solver_Volume_T {
	struct Solver_Volume_T volume; ///< The base \ref Solver_Volume_T.

	struct Multiarray_T* sol_coef_p; ///< The coefficients of the solution at a previous Runge-Kutta stage.
	struct Multiarray_T* rhs;        ///< The rhs terms.

	// Terms required for explicit runs.
	const struct const_Matrix_R* m_inv; ///< The inverse mass matrix.

	// Terms required for 2nd order PDE terms.
	struct Multiarray_T* grad_coef_v; ///< The volume contribution to the solution gradient coefficients.

	/// Linearization of \ref DG_Solver_Volume_T::grad_coef_v wrt \ref Solver_Volume_T::sol_coef.
	const struct const_Matrix_R* d_g_coef_v__d_s_coef[DIM];
};

/// \brief Constructor for a derived \ref DG_Solver_Volume_T.
void constructor_derived_DG_Solver_Volume_T
	(struct Volume* volume_ptr,   ///< Pointer to the volume.
	 const struct Simulation* sim ///< \ref Simulation.
	);

/// \brief Destructor for a derived \ref DG_Solver_Volume_T.
void destructor_derived_DG_Solver_Volume_T
	(struct Volume* volume_ptr ///< Pointer to the volume.
	);
