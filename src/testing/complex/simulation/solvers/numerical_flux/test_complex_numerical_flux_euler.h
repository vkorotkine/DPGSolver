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

#ifndef DPG__test_complex_numerical_flux_euler_h__INCLUDED
#define DPG__test_complex_numerical_flux_euler_h__INCLUDED
/** \file
 *  \brief Provides `complex` versions of functions defined in \ref numerical_flux_euler.h.
 */

struct Numerical_Flux_Input_c;
struct mutable_Numerical_Flux_c;

/// \brief `complex` version of \ref compute_Numerical_Flux_euler_lax_friedrichs.
void compute_Numerical_Flux_c_euler_lax_friedrichs
	(const struct Numerical_Flux_Input_c* num_flux_i, ///< See brief.
	 struct mutable_Numerical_Flux_c* num_flux        ///< See brief.
	);

/// \brief `complex` version of \ref compute_Numerical_Flux_euler_roe_pike.
void compute_Numerical_Flux_euler_roe_pike
	(const struct Numerical_Flux_Input_c* num_flux_i, ///< See brief.
	 struct mutable_Numerical_Flux_c* num_flux        ///< See brief.
	);

/// \brief `complex` version of \ref compute_Numerical_Flux_euler_roe_pike_jacobian.
void compute_Numerical_Flux_euler_roe_pike_jacobian
	(const struct Numerical_Flux_Input_c* num_flux_i, ///< See brief.
	 struct mutable_Numerical_Flux_c* num_flux        ///< See brief.
	);

#endif // DPG__test_complex_numerical_flux_euler_h__INCLUDED
