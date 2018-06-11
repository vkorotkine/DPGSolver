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
 *  \brief Provides templated functions relating to euler fluxes.
 */

#include "def_templates_flux.h"

struct Flux_Input_T;
struct mutable_Flux_T;

/** \brief Version of \ref compute_Flux_fptr_T computing the fluxes for the Euler equations.
 *
 *  The definitions of the Fluxes, Jacobians, and Hessians can be found in:
 *  - Fluxes:    \ref compute_Flux_euler_0;
 *  - Jacobians: \ref compute_Flux_euler_1;
 *  - Hessians:  \ref compute_Flux_euler_2.
 */
void compute_Flux_T_euler
	(const struct Flux_Input_T* flux_i, ///< See brief.
	 struct mutable_Flux_T* flux        ///< See brief.
	);

/** \brief Compute the square of the Velocity magnitude from the vector of uvw velocity components.
 *  \return See brief. */
Type compute_V2_from_uvw_T
	(const Type*const uvw ///< The array of velocity components.
	);

/** \brief Compute the square of the Velocity magnitude from the density and vector of rhouvw momentum components.
 *  \return See brief. */
Type compute_V2_from_rhouvw_T
	(const Type rho,         ///< The density.
	 const Type*const rhouvw ///< The array of momentum components.
	);

#include "undef_templates_flux.h"
