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

#ifndef DPG__flux_euler_h__INCLUDED
#define DPG__flux_euler_h__INCLUDED
/** \file
 *  \brief Provides functions relating to euler fluxes.
 */

struct Flux_Input;
struct mutable_Flux;

/// \brief `double` version of \ref compute_Flux_T_euler.
void compute_Flux_euler
	(const struct Flux_Input* flux_i, ///< \ref Flux_Input.
	 struct mutable_Flux* flux        ///< \ref Flux.
	);

/// \brief `double` version of \ref compute_Flux_T_euler_jacobian.
void compute_Flux_euler_jacobian
	(const struct Flux_Input* flux_i, ///< \ref Flux_Input.
	 struct mutable_Flux* flux        ///< \ref Flux.
	);

#endif // DPG__flux_euler_h__INCLUDED
