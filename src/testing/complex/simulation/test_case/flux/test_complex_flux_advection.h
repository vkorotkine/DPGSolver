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

#ifndef DPG__test_complex_flux_advection_h__INCLUDED
#define DPG__test_complex_flux_advection_h__INCLUDED
/** \file
 *  \brief Provides `complex` versions of containers and functions defined in \ref flux_advection.h.
 */

struct Flux_Input_c;
struct mutable_Flux_c;

/// \brief `complex` version of \ref compute_Flux_advection.
void compute_Flux_c_advection
	(const struct Flux_Input_c* flux_i, ///< \ref Flux_Input_c.
	 struct mutable_Flux_c* flux        ///< \ref Flux_c.
	);

/// \brief `complex` version of \ref compute_Flux_advection_jacobian.
void compute_Flux_c_advection_jacobian
	(const struct Flux_Input_c* flux_i, ///< \ref Flux_Input_c.
	 struct mutable_Flux_c* flux        ///< \ref Flux_c.
	);

#endif // DPG__test_complex_flux_advection_h__INCLUDED
