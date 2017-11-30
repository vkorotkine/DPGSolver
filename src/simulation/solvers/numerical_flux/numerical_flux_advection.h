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

#ifndef DPG__numerical_flux_advection_h__INCLUDED
#define DPG__numerical_flux_advection_h__INCLUDED
/** \file
 *  \brief Provides functions relating to linear advection numerical fluxes.
 */

struct Numerical_Flux_Input;
struct mutable_Numerical_Flux;

/// \brief `double` version of \ref compute_Numerical_Flux_T_advection_upwind.
void compute_Numerical_Flux_advection_upwind
	(const struct Numerical_Flux_Input* num_flux_i, ///< See brief.
	 struct mutable_Numerical_Flux* num_flux        ///< See brief.
	);

/// \brief `double` version of \ref compute_Numerical_Flux_T_advection_upwind_jacobian.
void compute_Numerical_Flux_advection_upwind_jacobian
	(const struct Numerical_Flux_Input* num_flux_i, ///< See brief.
	 struct mutable_Numerical_Flux* num_flux        ///< See brief.
	);

#endif // DPG__numerical_flux_advection_h__INCLUDED
