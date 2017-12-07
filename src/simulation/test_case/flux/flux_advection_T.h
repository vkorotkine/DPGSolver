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
 *  \brief Provides templated functions relating to linear advection fluxes.
 */

struct Flux_Input_T;
struct mutable_Flux_T;

/// \brief Version of \ref compute_Flux_fptr computing the fluxes for the linear advection equation.
void compute_Flux_T_advection
	(const struct Flux_Input_T* flux_i, ///< See brief.
	 struct mutable_Flux_T* flux        ///< See brief.
	);

/// \brief Version of \ref compute_Flux_fptr computing the fluxes and Jacobians for the linear advection equation.
void compute_Flux_T_advection_jacobian
	(const struct Flux_Input_T* flux_i, ///< See brief.
	 struct mutable_Flux_T* flux        ///< See brief.
	);
