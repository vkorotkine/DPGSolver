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

#ifndef DPG__test_complex_flux_euler_h__INCLUDED
#define DPG__test_complex_flux_euler_h__INCLUDED
/** \file
 *  \brief Provides `complex` versions of containers and functions defined in \ref flux_euler.h.
 */

struct Flux_Input_c;
struct mutable_Flux_;

/// \brief `complex` version of \ref compute_Flux_euler.
void compute_Flux_c_euler
	(const struct Flux_Input_c* flux_i, ///< See brief.
	 struct mutable_Flux_c* flux        ///< See brief.
	);

/// \brief `complex` version of \ref compute_Flux_euler_jacobian
void compute_Flux_c_euler_jacobian
	(const struct Flux_Input_c* flux_i, ///< See brief.
	 struct mutable_Flux_c* flux        ///< See brief.
	);

#endif // DPG__test_complex_flux_euler_h__INCLUDED
