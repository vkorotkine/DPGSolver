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
 *  \brief Provides templated functions relating to Navier-Stokes fluxes.
 */

struct Flux_Input_T;
struct mutable_Flux_T;

/** \brief Version of \ref compute_Flux_fptr_T computing the fluxes for the Navier-Stokes equations.
 *
 *  The implementation is based off of the fluxes as defined in Toro (Sections 1.3 and 1.4, \cite Toro2009).
 *
 *  \warning  Negated flux contributions are returned such that the identical treatment to inviscid flux contributions
 *            is applicable.
 *
 *  The definitions of the Fluxes, Jacobians, and Hessians can be found in:
 *  - Fluxes:    \ref compute_Flux_navier_stokes_0;
 *  - Jacobians: \ref compute_Flux_navier_stokes_1s and \ref compute_Flux_navier_stokes_1g.
 *  - Hessians:  \todo ref here.
 */
void compute_Flux_T_navier_stokes
	(const struct Flux_Input_T* flux_i, ///< See brief.
	 struct mutable_Flux_T* flux        ///< See brief.
	);
